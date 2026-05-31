"""
validation_full.jl — Julia values for all components and solvers.

Saves: results/julia_validation_full.csv
Compare with: python python/validation_full.py
Then diff: python python/compare_validation.py

Tests mirror python/validation_full.py exactly (same network, same parameters).
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Printf
using CSV       # stdlib-free fallback: write manually

# ── CSV writer ────────────────────────────────────────────────────────────────
rows = Vector{NamedTuple{(:test_id,:variable,:t,:value),
              Tuple{String,String,Int,Float64}}}()

function save(test_id::String, variable::String, t::Int, value::Real)
    push!(rows, (test_id=test_id, variable=variable,
                 t=t, value=round(Float64(value), digits=8)))
end

LP6 = [0.60, 0.57, 0.55, 0.54, 0.55, 0.60]

# ── Network builders ──────────────────────────────────────────────────────────
function base_3bus(; s_nom=1e9, x_l23=0.15)
    net = Network(baseMVA=100.0)
    add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
    add!(net, "Bus", "Bus2"; v_nom=380.0)
    add!(net, "Bus", "Bus3"; v_nom=380.0)
    add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,   s_nom=s_nom)
    add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,   s_nom=s_nom)
    add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=x_l23, s_nom=s_nom)
    return net
end

function add_gens_loads!(net; p_g1=400.0, p_g2=300.0, p_d2=200.0, p_d3=300.0)
    add!(net, "Generator", "G1"; bus="Bus1", p_nom=p_g1, marginal_cost=20.0)
    add!(net, "Generator", "G2"; bus="Bus2", p_nom=p_g2, marginal_cost=50.0)
    add!(net, "Load",      "D2"; bus="Bus2", p_set=p_d2)
    add!(net, "Load",      "D3"; bus="Bus3", p_set=p_d3)
end

# ─────────────────────────────────────────────────────────────────────────────
#  01  DC Power Flow
# ─────────────────────────────────────────────────────────────────────────────
print("01 DC Power Flow ... ")
net = base_3bus()
add!(net, "Generator", "G1"; bus="Bus1", p_nom=600.0, marginal_cost=20.0)
add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)

r = pf(net; method=:dc, verbose=false)
# θ is indexed by bus_names(net) order (sorted alphabetically)
bnames_sorted = bus_names(net)   # ["Bus1","Bus2","Bus3"]
for (i, b) in enumerate(bnames_sorted)
    save("01_dc_pf", "theta_$b", 0, r.θ[i])
end
# line_flows is a Vector of NamedTuple (name, from, to, P_MW, kind)
for lf in r.line_flows
    save("01_dc_pf", "flow_$(lf.name)", 0, lf.P_MW)
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  02  LOPF uncongested
# ─────────────────────────────────────────────────────────────────────────────
print("02 LOPF uncongested ... ")
net = base_3bus(); add_gens_loads!(net)
r = optimize(net; verbose=false)

save("02_lopf_free", "P_G1",      0, r.P_gen["G1"])
save("02_lopf_free", "P_G2",      0, r.P_gen["G2"])
save("02_lopf_free", "total_cost",0, r.total_cost)
for b in ["Bus1","Bus2","Bus3"]
    save("02_lopf_free", "LMP_$b", 0, r.lmp[b])
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  03  LOPF congested (s_nom=200 MW)
# ─────────────────────────────────────────────────────────────────────────────
print("03 LOPF congested ... ")
net = base_3bus(s_nom=200.0); add_gens_loads!(net)
r = optimize(net; verbose=false)

save("03_lopf_cong", "P_G1",      0, r.P_gen["G1"])
save("03_lopf_cong", "P_G2",      0, r.P_gen["G2"])
save("03_lopf_cong", "total_cost",0, r.total_cost)
for b in ["Bus1","Bus2","Bus3"]
    save("03_lopf_cong", "LMP_$b", 0, r.lmp[b])
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  04  LOPF with Transformer (tap_ratio=1)
# ─────────────────────────────────────────────────────────────────────────────
print("04 LOPF with Transformer ... ")
net = Network(baseMVA=100.0)
add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
add!(net, "Bus", "Bus2"; v_nom=380.0)
add!(net, "Bus", "Bus3"; v_nom=380.0)
add!(net, "Transformer", "T12";
     bus0="Bus1", bus1="Bus2",
     x=0.1, s_nom=500.0, tap_ratio=1.0, phase_shift=0.0)
add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e9)
add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e9)
add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
add!(net, "Generator", "G2"; bus="Bus2", p_nom=300.0, marginal_cost=50.0)
add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)
r = optimize(net; verbose=false)

save("04_lopf_trafo", "P_G1",      0, r.P_gen["G1"])
save("04_lopf_trafo", "P_G2",      0, r.P_gen["G2"])
save("04_lopf_trafo", "total_cost",0, r.total_cost)
# Transformer flow: DC approximation b_eff*(θ_from - θ_to)
bidx = bus_index(net)
tr   = net.transformers["T12"]
b_eff = net.baseMVA / (tr.x * tr.tap_ratio)
f_idx = bidx[tr.from_bus]; t_idx = bidx[tr.to_bus]
flow_T12 = b_eff * (r.θ[f_idx] - r.θ[t_idx])
save("04_lopf_trafo", "flow_T12", 0, flow_T12)
for b in ["Bus1","Bus2","Bus3"]
    save("04_lopf_trafo", "LMP_$b", 0, r.lmp[b])
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  05  LOPF with HVDC Link
# ─────────────────────────────────────────────────────────────────────────────
print("05 LOPF with HVDC Link ... ")
net = Network(baseMVA=100.0)
add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
add!(net, "Bus", "Bus2"; v_nom=380.0)
add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.2, s_nom=150.0)
add!(net, "Link", "HVDC"; bus0="Bus1", bus1="Bus2", p_nom=200.0, efficiency=0.97)
add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
add!(net, "Generator", "G2"; bus="Bus2", p_nom=200.0, marginal_cost=80.0)
add!(net, "Load", "D2"; bus="Bus2", p_set=350.0)
r = optimize(net; verbose=false)

save("05_lopf_link", "P_G1",      0, r.P_gen["G1"])
save("05_lopf_link", "P_G2",      0, r.P_gen["G2"])
save("05_lopf_link", "p_HVDC",    0, r.P_link["HVDC"])
save("05_lopf_link", "total_cost",0, r.total_cost)
for b in ["Bus1","Bus2"]
    save("05_lopf_link", "LMP_$b", 0, r.lmp[b])
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  06  LOPF with CO2 GlobalConstraint
# ─────────────────────────────────────────────────────────────────────────────
print("06 LOPF with CO2 cap ... ")
net = Network(baseMVA=100.0)
add!(net, "Bus",     "B1";   v_nom=380.0, slack=true)
add!(net, "Carrier", "coal"; co2_emissions=0.34)
add!(net, "Carrier", "gas";  co2_emissions=0.20)
add!(net, "Generator", "Coal"; bus="B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
add!(net, "Generator", "Gas";  bus="B1", p_nom=300.0, marginal_cost=50.0, carrier="gas")
add!(net, "Load",  "D1"; bus="B1", p_set=400.0)
add!(net, "GlobalConstraint", "co2_cap";
     constant=80.0, sense="<=",
     carrier_weightings=Dict("coal"=>0.34, "gas"=>0.20))
r = optimize(net; verbose=false)

save("06_lopf_co2", "P_Coal",    0, r.P_gen["Coal"])
save("06_lopf_co2", "P_Gas",     0, r.P_gen["Gas"])
save("06_lopf_co2", "total_cost",0, r.total_cost)
save("06_lopf_co2", "LMP_B1",   0, r.lmp["B1"])
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  07  Multi-period LOPF T=6, StorageUnit
# ─────────────────────────────────────────────────────────────────────────────
print("07 Multi-period LOPF T=6 + StorageUnit ... ")
net = base_3bus()
add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
add!(net, "Generator", "G2"; bus="Bus2", p_nom=300.0, marginal_cost=50.0)
add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)
add!(net, "StorageUnit", "Bat";
     bus="Bus2", p_nom=60.0, e_nom=240.0,
     efficiency_charge=0.9, efficiency_discharge=0.9,
     cyclic_state_of_charge=false, e_initial=120.0)

r = lopf_multiperiod(net; T=6, load_profile=LP6, wind_profile=LP6, verbose=false)

save("07_mp_storage", "total_cost", 0, r.total_cost)
for t in 1:6
    ti = t - 1   # 0-indexed to match Python iloc[t]
    save("07_mp_storage", "P_G1",    ti, r.gen_dispatch["G1"][t])
    save("07_mp_storage", "P_G2",    ti, r.gen_dispatch["G2"][t])
    save("07_mp_storage", "SoC_Bat", ti, r.soc["Bat"][t+1])   # t+1: Julia soc[1]=e_initial
    save("07_mp_storage", "p_ch",    ti, r.p_ch["Bat"][t])
    for b in ["Bus1","Bus2","Bus3"]
        save("07_mp_storage", "LMP_$b", ti, r.lmp["$b"][t])
    end
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  08  Unit Commitment T=6
# ─────────────────────────────────────────────────────────────────────────────
print("08 Unit Commitment T=6 ... ")
net = base_3bus()
add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
add!(net, "Generator", "G_peak";
     bus="Bus2", p_nom=300.0, marginal_cost=50.0,
     committable=true, min_up_time=2, min_down_time=1,
     startup_cost=1000.0, p_min_pu=0.3)
add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)

r = unit_commitment(net; T=6, load_profile=LP6, wind_profile=LP6,
                    compute_lmp=true, verbose=false)

save("08_uc", "total_cost", 0, r.total_cost)
for t in 1:6
    ti = t - 1
    save("08_uc", "P_G1",     ti, r.P_gen["G1"][t])
    save("08_uc", "P_G_peak", ti, r.P_gen["G_peak"][t])
    save("08_uc", "u_G_peak", ti, Float64(r.u["G_peak"][t]))
    for b in ["Bus1","Bus2","Bus3"]
        save("08_uc", "LMP_$b", ti, r.lmp["$b"][t])
    end
end
println("done")

# ─────────────────────────────────────────────────────────────────────────────
#  Write CSV (manual, no CSV.jl dependency)
# ─────────────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "..", "results", "julia_validation_full.csv")
open(out, "w") do f
    println(f, "test_id,variable,t,value")
    for r in rows
        println(f, "$(r.test_id),$(r.variable),$(r.t),$(r.value)")
    end
end
println("\n[OK] $(length(rows)) rows → $out")
println("Next: python python/compare_validation.py")
