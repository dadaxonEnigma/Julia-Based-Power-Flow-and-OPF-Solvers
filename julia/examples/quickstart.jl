# PowerFlowJulia — Quick Start
# Mirrors the PyPSA 3-bus tutorial exactly.
#
# PyPSA equivalent:
#   import pypsa
#   n = pypsa.Network()
#   n.add("Bus", "B1", v_nom=380); ...
#   n.optimize()

include("../src/PowerFlowJulia.jl")
using .PowerFlowJulia
using Printf

# ── Build network ──────────────────────────────────────────────────────────
net = Network(name="3-bus tutorial", baseMVA=100.0)

add!(net, "Bus", "Bus1"; v_nom=380.0, slack=true)
add!(net, "Bus", "Bus2"; v_nom=380.0)
add!(net, "Bus", "Bus3"; v_nom=380.0)

add!(net, "Line", "L1-2"; bus0="Bus1", bus1="Bus2", r=0.01, x=0.1, s_nom=200.0)
add!(net, "Line", "L1-3"; bus0="Bus1", bus1="Bus3", r=0.01, x=0.1, s_nom=200.0)
add!(net, "Line", "L2-3"; bus0="Bus2", bus1="Bus3", r=0.01, x=0.1, s_nom=200.0)

add!(net, "Generator", "G1"; bus="Bus1", p_nom=400.0, marginal_cost=20.0, carrier="gas")
add!(net, "Generator", "G2"; bus="Bus2", p_nom=300.0, marginal_cost=50.0, carrier="coal")

add!(net, "Load", "D2"; bus="Bus2", p_set=200.0)
add!(net, "Load", "D3"; bus="Bus3", p_set=300.0)

println(net)

# ── DC Power Flow ──────────────────────────────────────────────────────────
println("\n--- DC Power Flow ---")
pf_res = pf(net)                       # same as: dc_pf(net)

println("Bus angles (rad):")
for (b, θ) in zip(pf_res.buses, pf_res.θ)
    @printf("  %-10s  %.6f\n", b, θ)
end

# ── Linearized AC Power Flow ───────────────────────────────────────────────
println("\n--- Linearized AC Power Flow ---")
lac_res = pf(net, method=:lac)

println("Bus voltages (p.u.) and angles (rad):")
for (i, b) in enumerate(lac_res.buses)
    @printf("  %-10s  |V|=%.6f  θ=%.6f\n", b, lac_res.V_mag[i], lac_res.V_ang[i])
end

# ── Single-period LOPF ─────────────────────────────────────────────────────
println("\n--- Linear OPF (single period) ---")
opt_res = optimize(net)                # same as: lopf(net)

println("Generator dispatch:")
for (name, p) in opt_res.P_gen
    @printf("  %-10s  %.1f MW\n", name, p)
end
@printf("Total cost: %.0f €/h\n", opt_res.total_cost)

println("\nLocational Marginal Prices:")
for (bus, lmp) in opt_res.lmp
    @printf("  %-10s  %.2f €/MWh\n", bus, lmp)
end

# ── 24-hour Multi-period LOPF ─────────────────────────────────────────────
println("\n--- Multi-period LOPF (24h) ---")
add!(net, "StorageUnit", "Bat2";
     bus="Bus2", p_nom=80.0, e_nom=320.0,
     efficiency_charge=0.95, efficiency_discharge=0.95,
     cyclic_state_of_charge=true)

mp_res = optimize(net, T=24)           # same as: lopf_multiperiod(net, T=24)
@printf("24h total cost: %.0f €\n", mp_res.total_cost)
