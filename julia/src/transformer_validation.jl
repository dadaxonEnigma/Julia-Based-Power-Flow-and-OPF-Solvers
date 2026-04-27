"""
    transformer_validation.jl

Validates the transformer tap-ratio model against the IEEE 14-bus
MATPOWER reference solution.

Two DC PF runs are compared:
  (A) tap_ratio = 1.0 for all transformers  (ignoring tap)
  (B) tap_ratio from case14.m               (corrected)

Expected outcome: run B should show significantly smaller angle errors
at buses 6–14 which are downstream of the four off-nominal transformers.
"""

include("dispatch.jl")
using PowerModels
using Printf

CASE14 = joinpath(@__DIR__, "..", "..", "data", "case14.m")

# MATPOWER reference angles [degrees], buses 1–14
const REF_VA_DEG = [0.00, -4.98, -12.72, -10.33, -8.78,
                    -14.22, -13.37, -13.36, -14.94, -15.10,
                    -14.79, -15.07, -15.16, -16.04]

# ────────────────────────────────────────────────────────────────────────────
#  Parse case14.m with PowerModels.jl, build Network
# ────────────────────────────────────────────────────────────────────────────

function build_ieee14(; use_tap::Bool = true)
    data    = PowerModels.parse_file(CASE14)
    baseMVA = data["baseMVA"]
    net     = Network(name = use_tap ? "IEEE 14-bus (tap corrected)" :
                                       "IEEE 14-bus (tap ignored)",
                      baseMVA = baseMVA)

    # Buses
    for i in 1:14
        bus  = data["bus"]["$i"]
        add_bus!(net, "B$i",
                 v_nom    = bus["base_kv"],
                 v_mag_pu = bus["vm"],
                 v_ang    = bus["va"],
                 slack    = bus["bus_type"] == 3,
                 bus_type = bus["bus_type"])
    end

    # Branches: lines vs. transformers identified by tap ≠ 0 & ≠ 1
    branch_counter = Dict{String,Int}()
    for (_, br) in data["branch"]
        br["br_status"] == 0 && continue
        f = br["f_bus"]; t = br["t_bus"]

        # Unique name: "L f-t" or "T f-t", with counter for parallels
        is_trafo = get(br, "transformer", false) == true
        prefix   = is_trafo ? "T" : "L"
        base_key = "$prefix$f-$t"
        branch_counter[base_key] = get(branch_counter, base_key, 0) + 1
        cnt  = branch_counter[base_key]
        name = cnt == 1 ? base_key : "$(base_key)_$cnt"

        s_nom  = Inf   # thermal limits not used in DC PF validation

        if is_trafo
            tap = use_tap ? br["tap"] : 1.0
            add_transformer!(net, name, "B$f", "B$t",
                             r           = br["br_r"],
                             x           = br["br_x"],
                             s_nom       = s_nom,
                             tap_ratio   = tap,
                             phase_shift = br["shift"] * 180 / π)
        else
            add_line!(net, name, "B$f", "B$t",
                      r     = br["br_r"],
                      x     = br["br_x"],
                      b     = get(br, "b_fr", 0.0) + get(br, "b_to", 0.0),
                      s_nom = s_nom)
        end
    end

    # Generators: p_nom = actual dispatch pg (in MW) so dc_pf uses real injections.
    # p_max_pu = 1.0 → injection = p_nom * 1.0 = pg_MW.
    gen_counter = Dict{Int,Int}()
    for (_, gen) in data["gen"]
        gen["gen_status"] == 0 && continue
        b  = gen["gen_bus"]
        gen_counter[b] = get(gen_counter, b, 0) + 1
        cnt = gen_counter[b]
        pg_mw = gen["pg"] * baseMVA   # p.u. → MW
        add_generator!(net, "G$b" * (cnt > 1 ? "_$cnt" : ""), "B$b",
                       p_nom         = pg_mw,
                       p_max_pu      = 1.0,
                       marginal_cost = 0.0)
    end

    # Loads
    load_counter = Dict{Int,Int}()
    for (_, ld) in data["load"]
        b = ld["load_bus"]
        load_counter[b] = get(load_counter, b, 0) + 1
        cnt = load_counter[b]
        add_load!(net, "D$b" * (cnt > 1 ? "_$cnt" : ""), "B$b",
                  p_set = ld["pd"] * baseMVA,
                  q_set = ld["qd"] * baseMVA)
    end

    return net
end

# ────────────────────────────────────────────────────────────────────────────
#  Run comparison
# ────────────────────────────────────────────────────────────────────────────

net_tap    = build_ieee14(use_tap = true)
net_notap  = build_ieee14(use_tap = false)

println(net_tap)
println("\nTransformers in network:")
for (name, tr) in net_tap.transformers
    @printf("  %-8s  B%-2d → B%-2d  x=%.5f  tap=%.3f  shift=%.2f°\n",
            name, parse(Int, tr.from_bus[2:end]),
            parse(Int, tr.to_bus[2:end]),
            tr.x, tr.tap_ratio, tr.phase_shift)
end

res_tap   = dc_pf(net_tap,   verbose=false)
res_notap = dc_pf(net_notap, verbose=false)

# ── Comparison table ─────────────────────────────────────────────────────────
println()
println("="^78)
println("DC POWER FLOW — IEEE 14-bus: tap-corrected vs. tap-ignored")
println("="^78)
@printf("\n  %-5s  %10s  %10s  %10s  %10s  %10s\n",
        "Bus", "MATPOWER", "With tap", "|Δ| tap", "No tap", "|Δ| notap")
println("  " * "─"^68)

errs_tap   = Float64[]
errs_notap = Float64[]
for i in 1:14
    θ_tap   = res_tap.θ[res_tap.buses .== "B$i"][1]   * 180/π
    θ_notap = res_notap.θ[res_notap.buses .== "B$i"][1] * 180/π
    ref     = REF_VA_DEG[i]
    Δtap    = abs(θ_tap   - ref)
    Δnotap  = abs(θ_notap - ref)
    push!(errs_tap,   Δtap)
    push!(errs_notap, Δnotap)
    @printf("  B%-4d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n",
            i, ref, θ_tap, Δtap, θ_notap, Δnotap)
end
println("  " * "─"^68)
@printf("  %-5s  %10s  %10.4f  %10s  %10.4f  \n",
        "Max", "—", maximum(errs_tap), "—", maximum(errs_notap))
@printf("  %-5s  %10s  %10.4f  %10s  %10.4f  \n",
        "MAE", "—", sum(errs_tap)/14, "—", sum(errs_notap)/14)

improvement = (sum(errs_notap) - sum(errs_tap)) / sum(errs_notap) * 100
println()
@printf("Tap-ratio correction reduces MAE by %.1f%%.\n", improvement)
println()

# ── Transformer branch flows ──────────────────────────────────────────────────
println("Transformer branch flows (tap-corrected model):")
@printf("  %-8s  %12s  %12s\n", "Branch", "P from (MW)", "tap_ratio")
println("  " * "─"^36)
for f in res_tap.trafo_flows
    tr = net_tap.transformers[f.name]
    @printf("  %-8s  %12.2f  %12.3f\n", f.name, f.P_MW, tr.tap_ratio)
end
println()
println("="^78)
