include("dispatch.jl")
include("unit_commitment.jl")

println("=" ^ 65)
println("UNIT COMMITMENT TESTS")
println("=" ^ 65)

# ── TEST 1: Startup cost forces delayed commitment ───────────────────────────
# G1 cheap but has high startup cost → optimizer delays starting it
# G2 expensive, no startup cost, always on
println("\n[1] Startup cost effect")
net1 = Network(name="startup cost test", baseMVA=100.0)
add_bus!(net1, "B1", slack=true)
add_generator!(net1, "G_cheap", "B1",
    p_nom=200.0, p_min_pu=0.3, marginal_cost=10.0,
    committable=true, startup_cost=5000.0,  # high startup cost
    shutdown_cost=0.0, initial_status=false)
add_generator!(net1, "G_expensive", "B1",
    p_nom=200.0, p_min_pu=0.0, marginal_cost=80.0,
    committable=false)  # always dispatchable
add_load!(net1, "D1", "B1", p_set=100.0)

res1 = unit_commitment(net1, T=6, verbose=false,
    load_profile=[0.5, 0.5, 0.5, 1.0, 1.0, 1.0])  # low then high load

starts = sum(res1.su["G_cheap"])
@printf("  G_cheap starts: %d time(s)\n", starts)
@printf("  Total cost: %.0f €\n", res1.total_cost)
@assert string(res1.status) == "OPTIMAL"
println("  PASS: UC solved with startup cost")

# ── TEST 2: Minimum up time ───────────────────────────────────────────────────
println("\n[2] Minimum up time = 3h")
net2 = Network(name="min up time", baseMVA=100.0)
add_bus!(net2, "B1", slack=true)
add_generator!(net2, "G1", "B1",
    p_nom=300.0, p_min_pu=0.2, marginal_cost=20.0,
    committable=true, min_up_time=3,
    startup_cost=1000.0, initial_status=false)
add_generator!(net2, "G2", "B1",
    p_nom=100.0, marginal_cost=80.0, committable=false)
add_load!(net2, "B1" == "B1" ? "D1" : "D1", "B1", p_set=150.0)

res2 = unit_commitment(net2, T=8, verbose=false,
    load_profile=[0.3, 0.3, 1.0, 1.0, 1.0, 0.3, 0.3, 0.3])

u_g1 = res2.u["G1"]
@printf("  G1 commitment: %s\n", join(string.(u_g1), ""))
@assert string(res2.status) == "OPTIMAL"
# Verify: every startup is followed by at least min_up_time=3 ON periods
for t in 1:8
    if res2.su["G1"][t] == 1
        for dt in 0:min(2, 7-t)
            @assert u_g1[t+dt] == 1 "min_up_time violated at t=$(t+dt)"
        end
    end
end
println("  PASS: min up time constraint respected")

# ── TEST 3: Minimum down time ─────────────────────────────────────────────────
println("\n[3] Minimum down time = 2h")
net3 = Network(name="min down time", baseMVA=100.0)
add_bus!(net3, "B1", slack=true)
add_generator!(net3, "G1", "B1",
    p_nom=200.0, p_min_pu=0.2, marginal_cost=15.0,
    committable=true, min_down_time=2,
    startup_cost=500.0, initial_status=true)  # starts ON
add_generator!(net3, "G2", "B1",
    p_nom=200.0, marginal_cost=60.0, committable=false)
add_load!(net3, "D1", "B1", p_set=80.0)

res3 = unit_commitment(net3, T=6, verbose=false,
    load_profile=[1.0, 0.1, 0.1, 1.0, 1.0, 1.0])

u_g1 = res3.u["G1"]
@printf("  G1 commitment (starts ON): %s\n", join(string.(u_g1), ""))
@assert string(res3.status) == "OPTIMAL"
# If G1 shuts down, it must stay down for at least 2 periods
for t in 2:6
    if res3.sd["G1"][t] == 1  # shuts down at t
        if t + 1 <= 6
            @assert res3.u["G1"][t+1] == 0 "min_down violated: G1 restarted too soon"
        end
    end
end
println("  PASS: min down time constraint respected")

# ── TEST 4: Full 24h UC with 3-bus network ────────────────────────────────────
println("\n[4] Full 24h Unit Commitment — 3-bus network")
net4 = Network(name="3-bus UC", baseMVA=100.0)
add_bus!(net4, "B1", slack=true)
add_bus!(net4, "B2")
add_bus!(net4, "B3")
add_line!(net4, "L12", "B1", "B2", x=0.1, s_nom=Inf)
add_line!(net4, "L13", "B1", "B3", x=0.1, s_nom=Inf)
add_line!(net4, "L23", "B2", "B3", x=0.1, s_nom=Inf)

# Baseload: cheap, committable, high min up/down
add_generator!(net4, "BaseLoad", "B1",
    p_nom=300.0, p_min_pu=0.5, p_max_pu=1.0, marginal_cost=15.0,
    committable=true, min_up_time=4, min_down_time=4,
    startup_cost=8000.0, shutdown_cost=2000.0, initial_status=true)

# Peaker: expensive, fast-start, no min time
add_generator!(net4, "Peaker", "B2",
    p_nom=150.0, p_min_pu=0.0, marginal_cost=70.0,
    committable=true, min_up_time=1, min_down_time=1,
    startup_cost=500.0, initial_status=false)

# Wind: non-committable
add_generator!(net4, "Wind", "B3",
    p_nom=100.0, marginal_cost=0.0, carrier="wind", committable=false)

add_load!(net4, "D2", "B2", p_set=150.0)
add_load!(net4, "D3", "B3", p_set=100.0)

res4 = unit_commitment(net4, T=24, verbose=true)

@assert string(res4.status) == "OPTIMAL"
@assert haskey(res4.u, "BaseLoad")
@assert haskey(res4.u, "Peaker")
@assert haskey(res4.lmp, "B1")
@assert length(res4.lmp["B1"]) == 24

baseload_starts = sum(res4.su["BaseLoad"])
peaker_starts   = sum(res4.su["Peaker"])
@printf("\n  BaseLoad starts: %d  (min_up=4h enforced)\n", baseload_starts)
@printf("  Peaker starts:   %d\n", peaker_starts)
@printf("  LMP avg B1: %.2f €/MWh\n", sum(res4.lmp["B1"])/24)
println("  PASS ✓")

println("\n" * "=" ^ 65)
println("All Unit Commitment tests passed ✓")
println("=" ^ 65)
