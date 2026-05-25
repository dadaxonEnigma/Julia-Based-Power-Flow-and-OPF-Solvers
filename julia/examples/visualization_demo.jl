"""
Visualization demo — generates 5 plots from LOPF, multi-period LOPF, and UC results.
Saves PNGs to results/plots/.
"""

include(joinpath(@__DIR__, "..", "src", "PowerFlowJulia.jl"))
using .PowerFlowJulia

outdir = joinpath(@__DIR__, "..", "..", "results", "plots")
isdir(outdir) || mkpath(outdir)

LOAD_6 = [0.60, 0.57, 0.55, 0.54, 0.55, 0.60]

# ── 3-bus helper ─────────────────────────────────────────────────────────────
function make_3bus(; s_nom=Inf)
    net = Network(baseMVA=100.0)
    add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
    add!(net, "Bus", "Bus2"; v_nom=380.0)
    add!(net, "Bus", "Bus3"; v_nom=380.0)
    add!(net, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1, s_nom=s_nom)
    add!(net, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1, s_nom=s_nom)
    add!(net, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.1, s_nom=s_nom)
    add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0)
    add!(net, "Generator", "G2"; bus="Bus2", p_nom=300.0, marginal_cost=50.0)
    add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
    add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)
    return net
end

# ── 1. Network topology (congested case) ─────────────────────────────────────
println("1/5  Network topology ...")
net_con = make_3bus(s_nom=200.0)
r_con   = lopf(net_con; verbose=false)
p1 = plot_network(net_con; result=r_con,
                  title="3-Bus Network (s_nom=200 MW)",
                  savepath=joinpath(outdir, "network_topology.png"))
println("     saved network_topology.png")

# ── 2. Single-period LMP bar ──────────────────────────────────────────────────
println("2/5  LMP bar (congested LOPF) ...")
p2 = plot_lmp(r_con;
              title="LMP — Congested LOPF (s_nom=200 MW)",
              savepath=joinpath(outdir, "lmp_lopf_con.png"))
println("     saved lmp_lopf_con.png")

# ── 3. Multi-period dispatch stacked area ────────────────────────────────────
println("3/5  Multi-period dispatch ...")
net_mp = make_3bus()
add!(net_mp, "StorageUnit", "Bat";
     bus="Bus2", p_nom=60.0, e_nom=240.0,
     efficiency_charge=0.9, efficiency_discharge=0.9,
     cyclic_state_of_charge=false, e_initial=120.0)
r_mp = lopf_multiperiod(net_mp; T=6, load_profile=LOAD_6,
                         wind_profile=LOAD_6, verbose=false)
p3 = plot_dispatch(r_mp; net=net_mp,
                   title="Generator Dispatch — Multi-Period LOPF (T=6)",
                   savepath=joinpath(outdir, "dispatch_mp.png"))
println("     saved dispatch_mp.png")

# ── 4. Storage SoC ────────────────────────────────────────────────────────────
println("4/5  Storage SoC ...")
p4 = plot_soc(r_mp;
              title="Battery SoC — Multi-Period LOPF (T=6)",
              savepath=joinpath(outdir, "soc_mp.png"))
println("     saved soc_mp.png")

# ── 5. UC commitment schedule ─────────────────────────────────────────────────
println("5/5  UC schedule ...")
net_uc = Network(baseMVA=100.0)
add!(net_uc, "Bus", "Bus1"; v_nom=380.0, slack=true)
add!(net_uc, "Bus", "Bus2"; v_nom=380.0)
add!(net_uc, "Bus", "Bus3"; v_nom=380.0)
add!(net_uc, "Line", "L12"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1, s_nom=1e6)
add!(net_uc, "Line", "L13"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
add!(net_uc, "Line", "L23"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.1, s_nom=1e6)
add!(net_uc, "Generator", "G_base"; bus="Bus1", p_nom=300.0, marginal_cost=20.0)
add!(net_uc, "Generator", "G_peak1";
     bus="Bus2", p_nom=200.0, marginal_cost=55.0,
     committable=true, min_up_time=2, min_down_time=1,
     startup_cost=800.0, p_min_pu=0.3)
add!(net_uc, "Generator", "G_peak2";
     bus="Bus3", p_nom=150.0, marginal_cost=80.0,
     committable=true, min_up_time=1, min_down_time=1,
     startup_cost=500.0, p_min_pu=0.2)
add!(net_uc, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net_uc, "Load", "D3"; bus="Bus3", p_set=200.0)

LOAD_24 = PowerFlowJulia.DEFAULT_LOAD_PROFILE
r_uc = unit_commitment(net_uc; T=24, load_profile=LOAD_24,
                        wind_profile=LOAD_24, compute_lmp=true, verbose=false)
p5 = plot_uc_schedule(r_uc;
                      title="Unit Commitment Schedule (T=24h)",
                      savepath=joinpath(outdir, "uc_schedule.png"))
println("     saved uc_schedule.png")

println("\n[OK] All plots saved to results/plots/")
println("     network_topology.png  |  lmp_lopf_con.png  |  dispatch_mp.png")
println("     soc_mp.png            |  uc_schedule.png")
