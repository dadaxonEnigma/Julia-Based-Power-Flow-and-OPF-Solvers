using PowerModels
using Ipopt
using JuMP
using Printf


function network_to_powermodels(buses, lines, generators, loads, v_nom=380.0)
    
    baseMVA = 100.0
    
 
    data = Dict{String, Any}(
        "name" => "test_network",
        "baseMVA" => baseMVA,
        "per_unit" => true,          
        "source_type" => "unknown",
        "source_version" => "0",
        "bus" => Dict{String, Any}(),
        "load" => Dict{String, Any}(),
        "gen" => Dict{String, Any}(),
        "branch" => Dict{String, Any}(),
        "shunt" => Dict{String, Any}(),   
        "dcline" => Dict{String, Any}(),
        "storage" => Dict{String, Any}(),
        "switch" => Dict{String, Any}()  
    )
    
    for (i, bus_name) in enumerate(buses)
        data["bus"]["$i"] = Dict{String, Any}(
            "bus_i" => i,
            "bus_type" => 1,
            "vmax" => 1.1,
            "vmin" => 0.9,
            "va" => 0.0,
            "vm" => 1.0,
            "base_kv" => v_nom,
            "zone" => 1,
            "index" => i
        )
    end
    
   
    gen_idx = 1
    for (bus, (P, V_set)) in generators
        data["gen"]["$gen_idx"] = Dict{String, Any}(
            "gen_bus" => bus,
            "pg" => P / baseMVA,
            "qg" => 0.0,
            "qmax" => 1000.0 / baseMVA,
            "qmin" => -1000.0 / baseMVA,
            "pmax" => P * 2 / baseMVA,
            "pmin" => 0.0,
            "vg" => V_set,
            "mbase" => baseMVA,
            "gen_status" => 1,
            "index" => gen_idx
        )
        

        if gen_idx == 1
            data["bus"]["$bus"]["bus_type"] = 3
        end
        
        gen_idx += 1
    end
    
  
    load_idx = 1
    for (bus, (P, Q)) in loads
        data["load"]["$load_idx"] = Dict{String, Any}(
            "load_bus" => bus,
            "pd" => P / baseMVA,
            "qd" => Q / baseMVA,
            "status" => 1,
            "index" => load_idx
        )
        load_idx += 1
    end
    

    for (i, (from, to, r, x)) in enumerate(lines)
        z_base = v_nom^2 / baseMVA
        r_pu = r / z_base
        x_pu = x / z_base
        
        b_charge = 0.0  
        data["branch"]["$i"] = Dict{String, Any}(
            "f_bus" => from,
            "t_bus" => to,
            "br_r" => r_pu,
            "br_x" => x_pu,
            "br_b" => b_charge,
            "g_fr" => 0.0,             
            "b_fr" => b_charge / 2.0, 
            "g_to" => 0.0,             
            "b_to" => b_charge / 2.0,  
            "rate_a" => 1000.0,
            "rate_b" => 1000.0,
            "rate_c" => 1000.0,
            "tap" => 1.0,
            "shift" => 0.0,
            "br_status" => 1,
            "angmin" => -π/3,  
            "angmax" => π/3,  
            "index" => i
        )
    end
    
    return data
end


function solve_ac_pf_powermodels(buses, lines, generators, loads, line_names, v_nom=380.0)
    
    println("="^60)
    println("AC POWER FLOW - USING POWERMODELS.JL")
    println("="^60)
    
   
    println("\nStep 1: Converting network to PowerModels format...")
    pm_data = network_to_powermodels(buses, lines, generators, loads, v_nom)
    
    println("✓ Network converted")
    println("  Buses: $(length(pm_data["bus"]))")
    println("  Lines: $(length(pm_data["branch"]))")
    println("  Generators: $(length(pm_data["gen"]))")
    println("  Loads: $(length(pm_data["load"]))")
    
    
    println("\nStep 2: Solving AC Power Flow...")
    
    solver = optimizer_with_attributes(Ipopt.Optimizer, 
                                       "print_level" => 0,
                                       "sb" => "yes")
    
    result = solve_ac_pf(pm_data, solver)
    
    if result["termination_status"] == LOCALLY_SOLVED || 
       result["termination_status"] == OPTIMAL
        
        println("✓ AC Power Flow converged!")
        
        println("\n" * "="^60)
        println("RESULTS")
        println("="^60)
        
        sol = result["solution"]
        
        println("\n1. VOLTAGE MAGNITUDES (p.u.):")
        V_mag = Float64[]
        V_ang = Float64[]
        
        for i in 1:length(buses)
            vm = sol["bus"]["$i"]["vm"]
            va = sol["bus"]["$i"]["va"]
            push!(V_mag, vm)
            push!(V_ang, va)
            @printf("  %s: %.6f p.u.\n", buses[i], vm)
        end
        
        println("\n2. VOLTAGE ANGLES (radians):")
        for i in 1:length(buses)
            @printf("  %s: %.6f rad (%.6f°)\n", 
                    buses[i], V_ang[i], V_ang[i] * 180/π)
        end
        
        println("\n3. LINE FLOWS:")
        @printf("%-15s %15s %15s\n", "Line", "P (MW)", "Q (MVAr)")
        println("-"^50)

        line_P = Float64[]
        line_Q = Float64[]

        baseMVA = pm_data["baseMVA"]
        PowerModels.update_data!(pm_data, sol)
        branch_flows = PowerModels.calc_branch_flow_ac(pm_data)

        for (i, name) in enumerate(line_names)
            pf = branch_flows["branch"]["$i"]["pf"] * baseMVA   # p.u. → МВт
            qf = branch_flows["branch"]["$i"]["qf"] * baseMVA   # p.u. → МВАр
            push!(line_P, pf)
            push!(line_Q, qf)
            @printf("%-15s %15.6f %15.6f\n", name, pf, qf)
        end
        
        println("\n4. GENERATOR OUTPUT:")
        @printf("%-15s %15s %15s\n", "Generator", "P (MW)", "Q (MVAr)")
        println("-"^50)
        
        gen_P = Float64[]
        gen_Q = Float64[]
        
        for (idx, (bus, (P_nom, V))) in enumerate(generators)
            pg = sol["gen"]["$idx"]["pg"] * baseMVA   # p.u. → МВт
            qg = sol["gen"]["$idx"]["qg"] * baseMVA   # p.u. → МВАр
            push!(gen_P, pg)
            push!(gen_Q, qg)
            @printf("Gen %d (Bus %s) %15.6f %15.6f\n", idx, buses[bus], pg, qg)
        end
        
        println("\n" * "="^60)
        println("AC POWER FLOW COMPLETED!")
        println("="^60)
        
        return (
            V_mag = V_mag,
            V_ang = V_ang,
            line_P = line_P,
            line_Q = line_Q,
            gen_P = gen_P,
            gen_Q = gen_Q,
            converged = true
        )
        
    else
        println("⚠️  AC Power Flow did not converge!")
        println("Status: $(result["termination_status"])")
        return nothing
    end
end



println("\n\nCreating test network (same as PyPSA)...")

buses = ["Bus 0", "Bus 1", "Bus 2"]

lines = [
    (1, 2, 0.01, 0.1),
    (1, 3, 0.01, 0.1),
    (2, 3, 0.01, 0.1)
]
line_names = ["Line 0-1", "Line 0-2", "Line 1-2"]

generators = Dict(1 => (500.0, 1.0))

loads = Dict(
    2 => (300.0, 0.0),
    3 => (200.0, 0.0)
)

results = solve_ac_pf_powermodels(buses, lines, generators, loads, line_names, 380.0)

if results !== nothing
    println("\n\n" * "="^60)
    println("COMPARISON WITH PyPSA")
    println("="^60)
    
    println("\nExpected (PyPSA):")
    println("  V_mag: 1.0, 0.999982, 0.999984")
    println("  V_ang: 0.0, -0.000185, -0.000162")
    println("  Line P: 266.67, 233.34, -33.33")
    println("  Line Q: 0.048, 0.040, -0.002")
    
    println("\nActual (Julia/PowerModels):")
    @printf("  V_mag: %.6f, %.6f, %.6f\n", results.V_mag...)
    @printf("  V_ang: %.6f, %.6f, %.6f\n", results.V_ang...)
    @printf("  Line P: %.2f, %.2f, %.2f\n", results.line_P...)
    @printf("  Line Q: %.3f, %.3f, %.3f\n", results.line_Q...)
    
    println("\n\nACCURACY CHECK:")
    v_mag_diff = maximum(abs.(results.V_mag - [1.0, 0.999982, 0.999984]))
    v_ang_diff = maximum(abs.(results.V_ang - [0.0, -0.000185, -0.000162]))
    p_diff = maximum(abs.(results.line_P - [266.67, 233.34, -33.33]))
    
    @printf("Max V_mag difference: %.2e\n", v_mag_diff)
    @printf("Max V_ang difference: %.2e\n", v_ang_diff)
    @printf("Max P difference: %.2e MW\n", p_diff)
    
    if v_mag_diff < 1e-4 && v_ang_diff < 1e-5 && p_diff < 1.0
        println("\n✓✓✓ VALIDATION SUCCESSFUL! Results match PyPSA!")
    else
        println("\n⚠️  Some differences detected (but might be acceptable)")
    end
end