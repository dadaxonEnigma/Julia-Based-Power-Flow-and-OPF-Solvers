using PowerModels
using Ipopt
using JuMP
using HiGHS
using LinearAlgebra
using Printf
using Statistics

# ----------------------------------------------------------------
#  Paths
# ----------------------------------------------------------------
CASE14 = joinpath(@__DIR__, "..", "..", "data", "case14.m")

SILENT_IPOPT  = optimizer_with_attributes(Ipopt.Optimizer,
                    "print_level" => 0, "sb" => "yes")
SILENT_HIGHS  = optimizer_with_attributes(HiGHS.Optimizer,
                    "output_flag" => false)

# ----------------------------------------------------------------
#  Reference solution stored in case14.m (MATPOWER solved state)
#  Vm [p.u.] and Va [deg] for buses 1-14
# ----------------------------------------------------------------
const REF_VM = [1.060, 1.045, 1.010, 1.019, 1.020,
                1.070, 1.062, 1.090, 1.056, 1.051,
                1.057, 1.055, 1.050, 1.036]

const REF_VA_DEG = [0.00, -4.98, -12.72, -10.33, -8.78,
                   -14.22, -13.37, -13.36, -14.94, -15.10,
                   -14.79, -15.07, -15.16, -16.04]

# ----------------------------------------------------------------
#  Helper: print section header
# ----------------------------------------------------------------
header(s) = (println("\n" * "="^60); println(s); println("="^60))

# ================================================================
#  1. AC POWER FLOW  (PowerModels + Ipopt)
# ================================================================
header("IEEE 14-Bus: AC Power Flow  (PowerModels.jl + Ipopt)")

data_ac = PowerModels.parse_file(CASE14)
result_ac = solve_ac_pf(data_ac, SILENT_IPOPT)

if result_ac["termination_status"] in (LOCALLY_SOLVED, OPTIMAL)
    sol = result_ac["solution"]

    vm = [sol["bus"]["$i"]["vm"] for i in 1:14]
    va = [sol["bus"]["$i"]["va"] * 180/π for i in 1:14]

    println("\nBus  |  Vm Julia  Vm Ref  |ΔVm|    Va Julia  Va Ref  |ΔVa|")
    println("-"^65)
    for i in 1:14
        @printf("%3d  |  %.4f   %.4f  %.2e  | %8.4f  %8.4f  %.2e\n",
            i, vm[i], REF_VM[i], abs(vm[i]-REF_VM[i]),
            va[i], REF_VA_DEG[i], abs(va[i]-REF_VA_DEG[i]))
    end

    max_dvm = maximum(abs.(vm .- REF_VM))
    max_dva = maximum(abs.(va .- REF_VA_DEG))
    @printf("\nMax |ΔVm| = %.2e p.u.   Max |ΔVa| = %.2e deg\n", max_dvm, max_dva)
    println(max_dvm < 1e-2 && max_dva < 0.1 ?
        "✓  AC PF converged and matches MATPOWER reference." :
        "⚠  Larger than expected deviation from reference.")
else
    println("AC PF did not converge: $(result_ac["termination_status"])")
end

# ================================================================
#  2. DC POWER FLOW  (sparse linear solve)
# ================================================================
header("IEEE 14-Bus: DC Power Flow  (sparse linear solve)")

data_dc = PowerModels.parse_file(CASE14)
n = length(data_dc["bus"])

# Build susceptance matrix in per-unit
B = zeros(n, n)
for (_, br) in data_dc["branch"]
    br["br_status"] == 0 && continue
    i, j = br["f_bus"], br["t_bus"]
    b = 1.0 / br["br_x"]
    B[i,i] += b;  B[j,j] += b
    B[i,j] -= b;  B[j,i] -= b
end

# Net injection vector [p.u.]
P = zeros(n)
for (_, gen) in data_dc["gen"]
    gen["gen_status"] == 0 && continue
    P[gen["gen_bus"]] += gen["pg"]
end
for (_, ld) in data_dc["load"]
    P[ld["load_bus"]] -= ld["pd"]
end

# Solve reduced system (remove slack bus 1)
idx = 2:n
θ = zeros(n)
θ[idx] = B[idx, idx] \ P[idx]
θ_deg = θ .* (180/π)

println("\nBus  |  θ Julia (deg)  θ Ref (deg)   |Δθ|")
println("-"^50)
for i in 1:14
    @printf("%3d  |  %10.4f    %10.4f   %.2e\n",
        i, θ_deg[i], REF_VA_DEG[i], abs(θ_deg[i]-REF_VA_DEG[i]))
end

max_dθ = maximum(abs.(θ_deg .- REF_VA_DEG))
@printf("\nMax |Δθ| = %.4f deg  (DC vs full-AC reference)\n", max_dθ)
println("Note: IEEE 14-bus has 4 transformers (tap ≠ 1). DC PF ignores tap ratios → larger errors (~0.3–1.4°) at buses downstream of transformers (6–14).")

# ================================================================
#  3. LOPF  (JuMP + HiGHS)
# ================================================================
header("IEEE 14-Bus: LOPF  (JuMP + HiGHS)")

data_lp = PowerModels.parse_file(CASE14)

gen_buses  = sort([g["gen_bus"] for (_, g) in data_lp["gen"] if g["gen_status"]==1])
gen_pmax   = Dict(g["gen_bus"] => g["pmax"]
                  for (_, g) in data_lp["gen"] if g["gen_status"]==1)
gen_cost   = Dict(g["gen_bus"] => (length(g["cost"]) >= 2 ? g["cost"][end-1] : 1.0)
                  for (_, g) in data_lp["gen"] if g["gen_status"]==1)

P_load = zeros(n)
for (_, ld) in data_lp["load"]
    P_load[ld["load_bus"]] += ld["pd"]
end

# Rebuild B matrix
B2 = zeros(n, n)
sus_list = Tuple{Int,Int,Float64}[]
for (_, br) in data_lp["branch"]
    br["br_status"] == 0 && continue
    i, j = br["f_bus"], br["t_bus"]
    b = 1.0 / br["br_x"]
    push!(sus_list, (i, j, b))
    B2[i,i] += b;  B2[j,j] += b
    B2[i,j] -= b;  B2[j,i] -= b
end

model = Model(SILENT_HIGHS)

@variable(model, θ_v[1:n])
@variable(model, P_gen[bus in gen_buses],
          lower_bound = 0.0, upper_bound = gen_pmax[bus])

slack = minimum(gen_buses)
@constraint(model, θ_v[slack] == 0.0)

for k in 1:n
    P_inj = (k in gen_buses) ? P_gen[k] : 0.0
    @constraint(model,
        sum(B2[k,m]*θ_v[m] for m in 1:n) == P_inj - P_load[k])
end

@objective(model, Min,
    sum(gen_cost[bus] * P_gen[bus] for bus in gen_buses))

optimize!(model)

if termination_status(model) == OPTIMAL
    println("\nGenerator dispatch (LOPF):")
    @printf("%-8s  %10s  %10s  %10s\n", "Bus", "P (MW)", "Pmax (MW)", "Cost (€/MWh)")
    println("-"^45)
    # PowerModels scales cost by baseMVA: c_scaled = c_physical * baseMVA
    # So physical marginal cost = c_scaled / baseMVA
    baseMVA = 100.0
    let total_cost = 0.0
        for bus in gen_buses
            p          = value(P_gen[bus]) * baseMVA    # p.u. → MW
            pm         = gen_pmax[bus] * baseMVA
            c_physical = gen_cost[bus] / baseMVA        # descale to €/MWh
            total_cost += c_physical * p
            @printf("%-8d  %10.2f  %10.2f  %10.2f\n", bus, p, pm, c_physical)
        end
        @printf("\nTotal cost: %.2f €/h\n", total_cost)
    end
    @printf("Total load: %.2f MW\n", sum(P_load)*100)
else
    println("LOPF did not solve: $(termination_status(model))")
end

println("\n" * "="^60)
println("IEEE 14-Bus validation complete.")
println("="^60)
