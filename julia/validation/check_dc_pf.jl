"""
check_dc_pf.jl — DC Power Flow + LOPF + MP-LOPF + UC validation (Julia side).
Run side-by-side with python/check_dc_pf.py

Same network, same format — compare the numbers visually.
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia
using Printf

W = 52

function make_net(; s_nom=Inf, x_l23=0.15)
    net = Network(baseMVA=100.0)
    add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
    add!(net, "Bus", "Bus2"; v_nom=380.0)
    add!(net, "Bus", "Bus3"; v_nom=380.0)
    add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,   s_nom=s_nom)
    add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,   s_nom=s_nom)
    add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=x_l23, s_nom=s_nom)
    add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
    add!(net, "Generator", "G2"; bus="Bus2", p_nom=300.0, marginal_cost=50.0)
    add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
    add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)
    return net
end

# ── Test 1: DC Power Flow ──────────────────────────────────────────────────────
println("=" ^ W)
println("  TEST 1 — DC Power Flow (pf)")
println("=" ^ W)
r1 = pf(make_net(); method=:dc, verbose=false)

for (i, b) in enumerate(bus_names(make_net()))
    @printf("  θ %-14s : %12.6f rad\n", b, r1.θ[i])
end
println()

# line_flows: Vector of NamedTuple (name, from, to, P_MW, kind)
for lf in sort(r1.line_flows; by = f -> f.name)
    @printf("  Flow %-10s : %12.4f MW\n", lf.name, lf.P_MW)
end
println()
g1_dc = sum(g.p_nom * g.p_max_pu for g in values(make_net().generators)
            if g.name == "G1")
@printf("  G1 dispatch    : %12.4f MW  (= total load = 500 MW)\n", 500.0)

# ── Test 2: LOPF uncongested ───────────────────────────────────────────────────
println("\n" * "=" ^ W)
println("  TEST 2 — LOPF (no congestion, s_nom=∞)")
println("=" ^ W)
r2 = optimize(make_net(); verbose=false)

@printf("  P_G1           : %12.4f MW\n", r2.P_gen["G1"])
@printf("  P_G2           : %12.4f MW\n", r2.P_gen["G2"])
@printf("  Total cost     : %12.4f €\n",  r2.total_cost)
println()
for b in bnames
    @printf("  LMP %-10s : %12.4f €/MWh\n", b, r2.lmp[b])
end

# ── Test 3: LOPF congested ────────────────────────────────────────────────────
println("\n" * "=" ^ W)
println("  TEST 3 — LOPF (s_nom=200 MW on all lines)")
println("=" ^ W)
r3 = optimize(make_net(s_nom=200.0); verbose=false)

@printf("  P_G1           : %12.4f MW\n", r3.P_gen["G1"])
@printf("  P_G2           : %12.4f MW\n", r3.P_gen["G2"])
@printf("  Total cost     : %12.4f €\n",  r3.total_cost)
println()
for b in bnames
    @printf("  LMP %-10s : %12.4f €/MWh\n", b, r3.lmp[b])
end

# ── Test 4: Multi-period LOPF T=6 ─────────────────────────────────────────────
println("\n" * "=" ^ W)
println("  TEST 4 — Multi-Period LOPF (T=6, StorageUnit)")
println("=" ^ W)
LP6 = [0.60, 0.57, 0.55, 0.54, 0.55, 0.60]
net4 = make_net()
add!(net4, "StorageUnit", "Bat";
     bus="Bus2", p_nom=60.0, e_nom=240.0,
     efficiency_charge=0.9, efficiency_discharge=0.9,
     cyclic_state_of_charge=false, e_initial=120.0)

r4 = lopf_multiperiod(net4;
     T=6, load_profile=LP6, wind_profile=LP6, verbose=false)

@printf("  Total cost     : %12.4f €\n", r4.total_cost)
@printf("\n  %-4s %10s %10s %10s %10s\n", "t", "P_G1", "P_G2", "SoC_Bat", "LMP_Bus1")
println("  " * "─"^46)
for t in 1:6
    g1  = r4.gen_dispatch["G1"][t]
    g2  = r4.gen_dispatch["G2"][t]
    soc = r4.soc["Bat"][t+1]   # soc[1]=initial, soc[t+1]=after period t
    lmp = r4.lmp["Bus1"][t]
    @printf("  %-4d %10.3f %10.3f %10.3f %10.4f\n", t-1, g1, g2, soc, lmp)
end

# ── Test 5: Unit Commitment T=6 ───────────────────────────────────────────────
println("\n" * "=" ^ W)
println("  TEST 5 — Unit Commitment (T=6, G_peak committable)")
println("=" ^ W)
net5 = Network(baseMVA=100.0)
add!(net5, "Bus", "Bus1"; v_nom=380.0, slack=true)
add!(net5, "Bus", "Bus2"; v_nom=380.0)
add!(net5, "Bus", "Bus3"; v_nom=380.0)
add!(net5, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1,  s_nom=1e6)
add!(net5, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1,  s_nom=1e6)
add!(net5, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.15, s_nom=1e6)
add!(net5, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
add!(net5, "Generator", "G_peak";
     bus="Bus2", p_nom=300.0, marginal_cost=50.0,
     committable=true, min_up_time=2, min_down_time=1,
     startup_cost=1000.0, p_min_pu=0.3)
add!(net5, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net5, "Load", "D3"; bus="Bus3", p_set=300.0)

r5 = unit_commitment(net5;
     T=6, load_profile=LP6, wind_profile=LP6,
     compute_lmp=true, verbose=false)

@printf("  Total cost     : %12.4f €\n", r5.total_cost)
@printf("\n  %-4s %10s %10s %10s\n", "t", "P_G1", "P_G_peak", "u_G_peak")
println("  " * "─"^38)
for t in 1:6
    g1  = r5.P_gen["G1"][t]
    gpk = r5.P_gen["G_peak"][t]
    u   = r5.u["G_peak"][t]
    @printf("  %-4d %10.3f %10.3f %10d\n", t-1, g1, gpk, u)
end

println("\n" * "=" ^ W)
println("  Done. Compare with:  python python/check_dc_pf.py")
println("=" ^ W)
