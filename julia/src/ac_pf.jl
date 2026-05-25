"""
    ac_pf.jl

Full nonlinear AC Power Flow via PowerModels.jl + Ipopt.
Wraps the Network struct into PowerModels data format and calls solve_ac_pf.

This completes the pf() API parity with PyPSA's `network.pf()` (Newton-Raphson).
"""

@isdefined(Network) || include(joinpath(@__DIR__, "network.jl"))

using PowerModels
using Ipopt
import JuMP: optimizer_with_attributes

# ── Internal: convert Network → PowerModels data dict ────────────────────────
function _network_to_pm(net::Network)
    bnames = bus_names(net)
    bidx   = bus_index(net)
    n      = length(bnames)

    # Bus type: 3=slack, 2=PV (has generator), 1=PQ
    gen_buses  = Set(g.bus for g in values(net.generators) if g.carrier != "wind")
    slack_set  = Set(slack_buses(net))

    pm_buses = Dict{String,Any}()
    for (i, bname) in enumerate(bnames)
        bus_type = bname in slack_set ? 3 :
                   bname in gen_buses ? 2 : 1
        pm_buses["$i"] = Dict{String,Any}(
            "bus_i"    => i,
            "bus_type" => bus_type,
            "vm"       => net.buses[bname].v_mag_pu,
            "va"       => net.buses[bname].v_ang,
            "vmax"     => 1.1,
            "vmin"     => 0.9,
            "base_kv"  => net.buses[bname].v_nom,
            "zone"     => 1,
            "index"    => i,
        )
    end

    # Generators: set pg from p_nom*p_max_pu; large Q bounds
    pm_gens = Dict{String,Any}()
    for (gi, g) in enumerate(values(net.generators))
        g.carrier == "wind" && continue
        pm_gens["$gi"] = Dict{String,Any}(
            "gen_bus"    => bidx[g.bus],
            "pg"         => g.p_nom * g.p_max_pu / net.baseMVA,
            "qg"         => 0.0,
            "qmax"       =>  999.0,
            "qmin"       => -999.0,
            "pmax"       => g.p_nom * g.p_max_pu / net.baseMVA,
            "pmin"       => g.p_nom * g.p_min_pu / net.baseMVA,
            "vg"         => 1.0,
            "mbase"      => net.baseMVA,
            "gen_status" => 1,
            "cost"       => [0.0, g.marginal_cost],
            "ncost"      => 2,
            "index"      => gi,
        )
    end
    # Ensure at least one generator entry (PowerModels requires non-empty)
    isempty(pm_gens) && (pm_gens["1"] = Dict{String,Any}(
        "gen_bus" => 1, "pg" => 0.0, "qg" => 0.0,
        "qmax" => 0.0, "qmin" => 0.0, "pmax" => 0.0, "pmin" => 0.0,
        "vg" => 1.0, "mbase" => net.baseMVA, "gen_status" => 0,
        "cost" => [0.0, 0.0], "ncost" => 2, "index" => 1))

    # Loads
    pm_loads = Dict{String,Any}()
    for (li, l) in enumerate(values(net.loads))
        pm_loads["$li"] = Dict{String,Any}(
            "load_bus" => bidx[l.bus],
            "pd"       => l.p_set / net.baseMVA,
            "qd"       => l.q_set / net.baseMVA,
            "status"   => 1,
            "index"    => li,
        )
    end

    # Branches: Lines then Transformers
    pm_branches = Dict{String,Any}()
    bi = 1
    for l in values(net.lines)
        rate = isfinite(l.s_nom) ? l.s_nom / net.baseMVA : 9999.0
        pm_branches["$bi"] = Dict{String,Any}(
            "f_bus"     => bidx[l.from_bus],
            "t_bus"     => bidx[l.to_bus],
            "br_r"      => l.r,
            "br_x"      => l.x,
            "br_b"      => l.b,
            "g_fr"      => 0.0,
            "b_fr"      => l.b / 2.0,
            "g_to"      => 0.0,
            "b_to"      => l.b / 2.0,
            "rate_a"    => rate,
            "rate_b"    => rate,
            "rate_c"    => rate,
            "tap"       => 1.0,
            "shift"     => 0.0,
            "br_status" => 1,
            "angmin"    => -π / 3,
            "angmax"    =>  π / 3,
            "index"     => bi,
        )
        bi += 1
    end
    for tr in values(net.transformers)
        rate = isfinite(tr.s_nom) ? tr.s_nom / net.baseMVA : 9999.0
        pm_branches["$bi"] = Dict{String,Any}(
            "f_bus"     => bidx[tr.from_bus],
            "t_bus"     => bidx[tr.to_bus],
            "br_r"      => tr.r,
            "br_x"      => tr.x,
            "br_b"      => 0.0,
            "g_fr"      => 0.0,
            "b_fr"      => 0.0,
            "g_to"      => 0.0,
            "b_to"      => 0.0,
            "rate_a"    => rate,
            "rate_b"    => rate,
            "rate_c"    => rate,
            "tap"       => tr.tap_ratio,
            "shift"     => tr.phase_shift * π / 180,
            "br_status" => 1,
            "angmin"    => -π / 3,
            "angmax"    =>  π / 3,
            "index"     => bi,
        )
        bi += 1
    end

    return Dict{String,Any}(
        "name"           => net.name,
        "baseMVA"        => net.baseMVA,
        "per_unit"       => true,
        "source_type"    => "PowerFlowJulia",
        "source_version" => "0.1.0",
        "bus"            => pm_buses,
        "load"           => pm_loads,
        "gen"            => pm_gens,
        "branch"         => pm_branches,
        "shunt"          => Dict{String,Any}(),
        "dcline"         => Dict{String,Any}(),
        "storage"        => Dict{String,Any}(),
        "switch"         => Dict{String,Any}(),
    )
end

# ── Public solver ─────────────────────────────────────────────────────────────
"""
    ac_pf(net::Network; verbose=true) → NamedTuple

Full nonlinear AC Power Flow via PowerModels.jl + Ipopt.

Solves the full AC nodal power balance equations:
    P_i = Σⱼ |Vᵢ||Vⱼ|(Gᵢⱼ cos θᵢⱼ + Bᵢⱼ sin θᵢⱼ)
    Q_i = Σⱼ |Vᵢ||Vⱼ|(Gᵢⱼ sin θᵢⱼ − Bᵢⱼ cos θᵢⱼ)

Bus types (MATPOWER convention):
  3 = slack  (|V| and θ fixed, P and Q free)
  2 = PV     (|V| and P fixed, Q and θ free — buses with generators)
  1 = PQ     (P and Q fixed, |V| and θ free — load buses)

Returns:
  (V_mag, V_ang, P_flow, Q_flow, buses, converged)

PyPSA equivalent: `network.pf()` (Newton-Raphson)
"""
function ac_pf(net::Network; verbose=true)
    isempty(net.buses) && error("Network has no buses")
    isempty(slack_buses(net)) && error("Network has no slack bus defined")

    pm_data = _network_to_pm(net)
    bnames  = bus_names(net)
    n       = length(bnames)

    solver = optimizer_with_attributes(
        Ipopt.Optimizer,
        "print_level" => 0,
        "sb"          => "yes",
    )

    result    = PowerModels.solve_ac_pf(pm_data, solver)
    status    = result["termination_status"]
    converged = status in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)

    V_mag = fill(1.0, n)
    V_ang = zeros(n)
    if converged
        sol = result["solution"]
        for i in 1:n
            V_mag[i] = sol["bus"]["$i"]["vm"]
            V_ang[i] = sol["bus"]["$i"]["va"]
        end
    end

    # Branch flows
    line_list = collect(values(net.lines))
    P_flow    = Float64[]
    Q_flow    = Float64[]
    if converged && !isempty(line_list)
        PowerModels.update_data!(pm_data, result["solution"])
        flows = PowerModels.calc_branch_flow_ac(pm_data)
        for k in 1:length(line_list)
            push!(P_flow, get(get(flows["branch"], "$k", Dict()), "pf", 0.0) * net.baseMVA)
            push!(Q_flow, get(get(flows["branch"], "$k", Dict()), "qf", 0.0) * net.baseMVA)
        end
    end

    if verbose
        println("="^60)
        println("AC POWER FLOW (PowerModels.jl + Ipopt) — $(net.name)")
        println(converged ? "  Status: CONVERGED" : "  Status: DID NOT CONVERGE ($status)")
        println("="^60)
        if converged
            @printf("\n  %-14s  %10s  %10s\n", "Bus", "|V| (p.u.)", "θ (deg)")
            println("  " * "─"^38)
            for (i, b) in enumerate(bnames)
                @printf("  %-14s  %10.6f  %10.4f\n", b, V_mag[i], V_ang[i] * 180 / π)
            end
            if !isempty(P_flow)
                println("\n  Line flows:")
                @printf("  %-14s  %10s  %10s\n", "Line", "P (MW)", "Q (MVAr)")
                println("  " * "─"^38)
                for (k, l) in enumerate(line_list)
                    @printf("  %-14s  %10.4f  %10.4f\n", l.name, P_flow[k], Q_flow[k])
                end
            end
        end
        println("="^60)
    end

    return (V_mag=V_mag, V_ang=V_ang, P_flow=P_flow, Q_flow=Q_flow,
            buses=bnames, converged=converged)
end
