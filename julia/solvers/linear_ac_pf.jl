using LinearAlgebra
using Printf

"""
    solve_linear_ac_pf(buses, lines, generators, loads; baseMVA, v_nom)

Linearized AC Power Flow (LACPF).

Unlike DC PF (which ignores resistance and reactive power), LACPF linearises
the full AC power flow equations around the flat start (|V|=1 p.u., θ=0 rad)
using a first-order Taylor expansion. It solves simultaneously for:
  - voltage angle deviations  Δθ  [rad]
  - voltage magnitude deviations  Δ|V| [p.u.]

The linearised system (excluding slack bus) is:
    [ B'  -G'] [Δθ ]   [P_inj / baseMVA]
    [-G'  -B'] [Δ|V|] = [Q_inj / baseMVA]

where G' + jB' is the reduced nodal admittance matrix.

# Arguments
- `buses`      : vector of bus names
- `lines`      : vector of (from, to, r_pu, x_pu) — r,x in per-unit
- `generators` : Dict(bus => (P_MW, Q_MVAr))  (Q=0 if not reactive-controlled)
- `loads`      : Dict(bus => (P_MW, Q_MVAr))
- `baseMVA`    : system base (default 100 MVA)
- `v_nom`      : nominal voltage kV (informational only; r,x already p.u.)
"""
function solve_linear_ac_pf(buses, lines, generators, loads;
                              baseMVA = 100.0,
                              v_nom   = 380.0,
                              verbose = true)

    n = length(buses)

    # ── Build admittance matrix Y = G + jB ──────────────────────
    # r, x are in physical Ω → convert to p.u. via z_base = v_nom²/baseMVA
    z_base = v_nom^2 / baseMVA          # e.g. 380²/100 = 1444 Ω

    G = zeros(n, n)
    B = zeros(n, n)

    for (from, to, r, x) in lines
        r_pu  = r / z_base               # per-unit resistance
        x_pu  = x / z_base               # per-unit reactance
        denom = r_pu^2 + x_pu^2         # |z_pu|²

        g_ij  =  r_pu / denom            # conductance [p.u.]
        b_ij  =  x_pu / denom            # |susceptance| [p.u.] (positive)

        # Diagonal (shunt): sum of incident branch admittances
        G[from, from] += g_ij;  G[to, to]   += g_ij
        B[from, from] += b_ij;  B[to, to]   += b_ij   # B_ii > 0

        # Off-diagonal: negative admittance
        G[from, to]   -= g_ij;  G[to, from] -= g_ij
        B[from, to]   -= b_ij;  B[to, from] -= b_ij   # B_ij < 0
    end

    # Scale to MW units: B [MW/rad], G [MW/p.u.-V]
    G_mw = G .* baseMVA
    B_mw = B .* baseMVA

    # ── Net power injections [p.u.] ──────────────────────────────
    P = zeros(n)   # active   [MW] → will be divided by baseMVA
    Q = zeros(n)   # reactive [MVAr]

    for (bus, (p, q)) in generators
        P[bus] += p
        Q[bus] += q
    end
    for (bus, (p, q)) in loads
        P[bus] -= p
        Q[bus] -= q
    end

    # ── Remove slack bus (bus 1) → reduce to (n-1) × (n-1) ─────
    idx = 2:n    # non-slack indices

    B_r = B_mw[idx, idx]    # reduced B [MW/rad]
    G_r = G_mw[idx, idx]    # reduced G [MW/p.u.-V]
    P_r = P[idx]             # [MW]
    Q_r = Q[idx]             # [MVAr]

    # ── Assemble 2(n-1) × 2(n-1) system ─────────────────────────
    #   [ B_r  -G_r] [Δθ ]   [P_r]
    #   [-G_r  -B_r] [Δ|V|] = [Q_r]
    nm = n - 1
    A = [ B_r  -G_r ;
         -G_r  -B_r ]

    rhs = [P_r; Q_r]

    # ── Solve ────────────────────────────────────────────────────
    x_sol = A \ rhs

    Δθ_red = x_sol[1:nm]        # [rad]
    Δv_red = x_sol[nm+1:end]    # [p.u.]

    # ── Reconstruct full vectors ─────────────────────────────────
    θ   = zeros(n);      θ[idx]   = Δθ_red
    Δv  = zeros(n);      Δv[idx]  = Δv_red
    V_mag = 1.0 .+ Δv           # |V| = 1 + Δ|V|  [p.u.]

    # ── Line flows ───────────────────────────────────────────────
    P_flow = Float64[]
    Q_flow = Float64[]

    for (from, to, r, x) in lines
        r_pu  = r / z_base
        x_pu  = x / z_base
        denom = r_pu^2 + x_pu^2
        g_ij  = r_pu / denom     # conductance p.u.
        b_ij  = x_pu / denom     # |susceptance| p.u.

        Δθ_ft = θ[from] - θ[to]
        ΔV_ft = Δv[from] - Δv[to]

        # Linearised branch flows (1st order, flat start):
        #   P_ft ≈ (b_ij * Δθ_ft + g_ij * ΔV_ft) * baseMVA
        #   Q_ft ≈ (b_ij * ΔV_ft - g_ij * Δθ_ft) * baseMVA
        p_ft = ( b_ij * Δθ_ft + g_ij * ΔV_ft) * baseMVA
        q_ft = ( b_ij * ΔV_ft - g_ij * Δθ_ft) * baseMVA

        push!(P_flow, p_ft)
        push!(Q_flow, q_ft)
    end

    # ── Print results ────────────────────────────────────────────
    if verbose
        println("="^60)
        println("LINEARIZED AC POWER FLOW  (Julia)")
        println("="^60)

        println("\n1. VOLTAGE MAGNITUDES AND ANGLES:")
        @printf("%-10s %12s %12s %12s\n",
                "Bus", "|V| (p.u.)", "θ (rad)", "θ (deg)")
        println("-"^50)
        for i in 1:n
            @printf("%-10s %12.6f %12.6f %12.4f\n",
                    buses[i], V_mag[i], θ[i], θ[i]*180/π)
        end

        println("\n2. LINE FLOWS:")
        @printf("%-15s %12s %12s\n", "Line", "P (MW)", "Q (MVAr)")
        println("-"^42)
        for (k, (from, to, r, x)) in enumerate(lines)
            @printf("Bus%d → Bus%d    %12.4f %12.4f\n",
                    from, to, P_flow[k], Q_flow[k])
        end
    end

    return (
        V_mag   = V_mag,
        V_ang   = θ,
        P_flow  = P_flow,
        Q_flow  = Q_flow
    )
end

# ================================================================
#  Main: 3-bus validation against DC PF and full AC PF reference
# ================================================================
if abspath(PROGRAM_FILE) == @__FILE__

println("\n" * "="^60)
println("3-BUS VALIDATION: DC PF vs LINEARIZED AC PF vs Full AC PF")
println("="^60)

buses = ["Bus 0", "Bus 1", "Bus 2"]
lines = [(1, 2, 0.01, 0.1),
         (1, 3, 0.01, 0.1),
         (2, 3, 0.01, 0.1)]

generators = Dict(1 => (500.0, 0.0))
loads      = Dict(2 => (300.0, 0.0),
                  3 => (200.0, 0.0))

res = solve_linear_ac_pf(buses, lines, generators, loads)

println("\n" * "="^60)
println("COMPARISON: Linearized AC PF vs Reference (PyPSA full AC PF)")
println("="^60)

# PyPSA full AC PF reference values (from validation in Chapter 3)
ref_vm  = [1.000000, 0.999982, 0.999984]
ref_va  = [0.000000, -0.000185, -0.000162]
ref_pf  = [266.67, 233.34, -33.33]
ref_qf  = [0.048, 0.040, -0.002]

println("\nVoltage magnitudes [p.u.]:")
@printf("%-8s %12s %12s %12s\n", "Bus", "LACPF", "Full AC", "|Δ|")
for i in 1:3
    @printf("%-8s %12.6f %12.6f %12.2e\n",
            buses[i], res.V_mag[i], ref_vm[i],
            abs(res.V_mag[i] - ref_vm[i]))
end

println("\nVoltage angles [rad]:")
@printf("%-8s %12s %12s %12s\n", "Bus", "LACPF", "Full AC", "|Δ|")
for i in 1:3
    @printf("%-8s %12.6f %12.6f %12.2e\n",
            buses[i], res.V_ang[i], ref_va[i],
            abs(res.V_ang[i] - ref_va[i]))
end

println("\nLine active power flows [MW]:")
line_names = ["L 0-1", "L 0-2", "L 1-2"]
# DC PF reference (from dc_power_flow.jl validation, baseMVA=100 convention)
dc_pf_flows = [266.67, 233.33, -33.33]  # same as AC for this lossless-like case
@printf("%-8s %10s %10s %10s %10s\n", "Line", "LACPF P", "Full AC P", "|ΔP|", "DC PF P")
for i in 1:3
    @printf("%-8s %10.4f %10.4f %10.4f %10.4f\n",
            line_names[i], res.P_flow[i], ref_pf[i],
            abs(res.P_flow[i] - ref_pf[i]), dc_pf_flows[i])
end

println("\nLine reactive power flows [MVAr]:")
@printf("%-8s %10s %10s %10s\n", "Line", "LACPF Q", "Full AC Q", "|ΔQ|")
for i in 1:3
    @printf("%-8s %10.4f %10.4f %10.4f\n",
            line_names[i], res.Q_flow[i], ref_qf[i],
            abs(res.Q_flow[i] - ref_qf[i]))
end

max_dvm = maximum(abs.(res.V_mag .- ref_vm))
max_dva = maximum(abs.(res.V_ang .- ref_va))
max_dp  = maximum(abs.(res.P_flow .- ref_pf))

println("\n" * "─"^50)
@printf("Max |ΔVm|:  %.2e p.u.\n", max_dvm)
@printf("Max |Δθ|:   %.2e rad\n",  max_dva)
@printf("Max |ΔP|:   %.4f MW\n",   max_dp)
println("─"^50)
println("Note: DC PF ignores R → no reactive power, Vm=1 always.")
println("LACPF recovers approximate Q flows and voltage deviations.")
println("="^60)

end # if abspath
