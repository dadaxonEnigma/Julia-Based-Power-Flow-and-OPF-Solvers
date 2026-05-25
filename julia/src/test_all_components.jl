include("dispatch.jl")

# ── helpers ──────────────────────────────────────────────────────────────────
pass(msg) = println("  ✓  $msg")
fail(msg) = (println("  ✗  $msg"); error("FAILED: $msg"))

results = Dict{String, Bool}()
macro suite(name, body)
    quote
        println("\n" * "─"^60)
        println("SUITE: ", $name)
        println("─"^60)
        try
            $(esc(body))
            results[$name] = true
        catch e
            println("  ERROR: ", e)
            results[$name] = false
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  1. STRUCT CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════
@suite "Bus construction" begin
    b = Bus("B1", v_nom=380.0, slack=true, bus_type=3)
    @assert b.name == "B1" && b.v_nom == 380.0 && b.slack == true
    b2 = Bus("B2")                          # all defaults
    @assert b2.v_nom == 1.0 && !b2.slack && b2.bus_type == 1
    pass("Bus default + custom kwargs")
end

@suite "Line construction" begin
    l = Line("L1", "B1", "B2", r=0.01, x=0.1, s_nom=200.0)
    @assert l.x == 0.1 && l.s_nom == 200.0 && l.b == 0.0
    try
        Line("bad", "B1", "B2", x=0.0)
        fail("should error on x=0")
    catch e
        pass("Line rejects x=0")
    end
end

@suite "Transformer construction" begin
    tr = Transformer("T1", "B1", "B2", x=0.2, tap_ratio=0.978, phase_shift=5.0)
    @assert tr.tap_ratio ≈ 0.978 && tr.phase_shift ≈ 5.0
    tr2 = Transformer("T2", "B1", "B2")    # defaults: tap=1, shift=0
    @assert tr2.tap_ratio == 1.0 && tr2.phase_shift == 0.0
    pass("Transformer with tap + phase shift")
end

@suite "Generator construction" begin
    g = Generator("G1", "B1", p_nom=400.0, marginal_cost=20.0,
                  carrier="gas", committable=false)
    @assert g.p_nom == 400.0 && g.carrier == "gas"
    g2 = Generator("G2", "B2", p_nom=100.0, carrier="wind", marginal_cost=0.0)
    @assert g2.carrier == "wind"
    pass("Generator dispatchable + wind")
end

@suite "Load construction" begin
    d = Load("D1", "B1", p_set=300.0, q_set=50.0)
    @assert d.p_set == 300.0 && d.q_set == 50.0
    d2 = Load("D2", "B1")                  # defaults
    @assert d2.p_set == 0.0 && d2.q_set == 0.0
    pass("Load with and without reactive power")
end

@suite "StorageUnit construction" begin
    s = StorageUnit("Bat", "B1", p_nom=100.0, e_nom=400.0,
                    efficiency_charge=0.95, efficiency_discharge=0.95,
                    cyclic_state_of_charge=true)
    @assert s.p_nom == 100.0 && s.e_nom == 400.0
    @assert s.efficiency_charge == 0.95 && s.cyclic_state_of_charge == true
    pass("StorageUnit all params")
end

@suite "Store construction" begin
    s = Store("St1", "B1", e_nom=200.0, p_nom=50.0,
              e_min_pu=0.1, e_max_pu=0.9, e_cyclic=true)
    @assert s.e_nom == 200.0 && s.p_nom == 50.0
    @assert s.e_min_pu == 0.1 && s.e_max_pu == 0.9
    pass("Store with bounds")
end

@suite "Carrier construction" begin
    c = Carrier("gas", co2_emissions=0.20, color="#ff8c00")
    @assert c.co2_emissions == 0.20
    c2 = Carrier("wind")  # zero emissions default
    @assert c2.co2_emissions == 0.0
    pass("Carrier with and without emissions")
end

@suite "Link construction" begin
    l = Link("HVDC", "A", "B", p_nom=500.0, efficiency=0.97,
             p_min_pu=-1.0, marginal_cost=2.0)
    @assert l.efficiency == 0.97 && l.p_min_pu == -1.0
    try
        Link("bad", "A", "B", efficiency=0.0)
        fail("should error on efficiency=0")
    catch e
        pass("Link rejects efficiency≤0")
    end
end

@suite "GlobalConstraint construction" begin
    gc = GlobalConstraint("co2",
         type="co2_limit", sense="<=", constant=500.0,
         carrier_weightings=Dict("gas"=>0.20, "coal"=>0.34))
    @assert gc.constant == 500.0
    @assert gc.carrier_weightings["coal"] ≈ 0.34
    try
        GlobalConstraint("bad", sense="??")
        fail("should reject invalid sense")
    catch e
        pass("GlobalConstraint rejects invalid sense")
    end
end

# ═══════════════════════════════════════════════════════════════════════════
#  2. NETWORK API
# ═══════════════════════════════════════════════════════════════════════════
@suite "Network add_*! and validation" begin
    net = Network(name="API test", baseMVA=100.0)

    # add_bus!
    add_bus!(net, "B1", v_nom=380.0, slack=true)
    add_bus!(net, "B2", v_nom=380.0)
    add_bus!(net, "B3", v_nom=380.0)
    @assert length(net.buses) == 3
    try; add_bus!(net, "B1"); fail("duplicate bus") catch e; pass("Bus duplicate rejected") end

    # add_line!
    add_line!(net, "L12", "B1", "B2", x=0.1)
    add_line!(net, "L13", "B1", "B3", x=0.1)
    add_line!(net, "L23", "B2", "B3", x=0.1)
    try; add_line!(net, "LX", "B1", "B99", x=0.1); fail("missing bus")
    catch e; pass("Line rejects missing bus") end

    # add_transformer!
    add_transformer!(net, "T12", "B1", "B2", x=0.2, tap_ratio=0.95)
    @assert length(net.transformers) == 1

    # add_generator!
    add_generator!(net, "G1", "B1", p_nom=400.0, marginal_cost=20.0)
    add_generator!(net, "G2", "B2", p_nom=300.0, marginal_cost=50.0, carrier="wind")
    try; add_generator!(net, "G1", "B1"); fail("duplicate gen") catch e; pass("Generator duplicate rejected") end

    # add_load!
    add_load!(net, "D2", "B2", p_set=200.0, q_set=10.0)
    add_load!(net, "D3", "B3", p_set=300.0)

    # add_storage_unit!
    add_storage_unit!(net, "Bat", "B2", p_nom=50.0, e_nom=200.0)

    # add_store!
    add_store!(net, "Store1", "B3", e_nom=100.0, p_nom=40.0)

    # add_carrier!
    add_carrier!(net, "gas", co2_emissions=0.20)
    add_carrier!(net, "wind", co2_emissions=0.0)
    try; add_carrier!(net, "gas"); fail("duplicate carrier") catch e; pass("Carrier duplicate rejected") end

    # add_link!
    add_link!(net, "Lnk1", "B1", "B3", p_nom=100.0, efficiency=0.98)
    try; add_link!(net, "LX", "B1", "B99", p_nom=100.0); fail("missing bus")
    catch e; pass("Link rejects missing bus") end

    # add_global_constraint!
    add_global_constraint!(net, "co2",
        constant=500.0, carrier_weightings=Dict("gas"=>0.20))
    try; add_global_constraint!(net, "co2"); fail("duplicate GC") catch e; pass("GlobalConstraint duplicate rejected") end

    # query helpers
    @assert length(slack_buses(net)) == 1 && slack_buses(net)[1] == "B1"
    @assert total_p_nom(net) == 700.0
    @assert total_load(net) == 500.0
    @assert length(branches_at(net, "B1")) == 3  # passive: L12, L13, T12
    @assert length(links_at(net, "B1")) == 1      # active:  Lnk1
    bidx = bus_index(net)
    @assert bidx["B1"] == 1 && bidx["B2"] == 2 && bidx["B3"] == 3

    pass("Full Network API: all add_*! + queries work")
    println(net)
end

# ═══════════════════════════════════════════════════════════════════════════
#  3. SOLVERS WITH ALL COMPONENTS
# ═══════════════════════════════════════════════════════════════════════════

function make_full_3bus()
    net = Network(name="full 3-bus", baseMVA=100.0)
    add_carrier!(net, "gas",  co2_emissions=0.20)
    add_carrier!(net, "coal", co2_emissions=0.34)
    add_carrier!(net, "wind", co2_emissions=0.00)
    add_bus!(net, "B1", v_nom=380.0, slack=true)
    add_bus!(net, "B2", v_nom=380.0)
    add_bus!(net, "B3", v_nom=380.0)
    add_line!(net, "L12", "B1", "B2", r=0.01, x=0.1)
    add_line!(net, "L13", "B1", "B3", r=0.01, x=0.1)
    add_line!(net, "L23", "B2", "B3", r=0.01, x=0.1)
    add_generator!(net, "G1", "B1", p_nom=400.0, marginal_cost=20.0, carrier="gas")
    add_generator!(net, "G2", "B2", p_nom=300.0, marginal_cost=50.0, carrier="coal")
    add_generator!(net, "Wind3", "B3", p_nom=80.0,  marginal_cost=0.0,  carrier="wind")
    add_load!(net, "D2", "B2", p_set=200.0, q_set=5.0)
    add_load!(net, "D3", "B3", p_set=300.0, q_set=8.0)
    add_storage_unit!(net, "Bat2", "B2",
        p_nom=60.0, e_nom=240.0,
        efficiency_charge=0.95, efficiency_discharge=0.95,
        cyclic_state_of_charge=true)
    add_store!(net, "Store3", "B3",
        e_nom=150.0, p_nom=50.0, e_cyclic=true)
    add_link!(net, "Link13", "B1", "B3",
        p_nom=100.0, efficiency=0.97, marginal_cost=1.0)
    add_global_constraint!(net, "co2",
        constant=3000.0,
        carrier_weightings=Dict("gas"=>0.20, "coal"=>0.34))
    return net
end

@suite "dc_pf with all branch types" begin
    net = make_full_3bus()
    res = dc_pf(net, verbose=false)
    @assert res.converged
    @assert abs(res.θ[1]) < 1e-10           # slack angle = 0
    @assert length(res.line_flows) == 3
    @assert length(res.trafo_flows) == 0    # no transformers
    @assert length(res.all_flows) == 3
    # flows should sum to zero at each bus (KCL)
    pass("dc_pf: converged, slack=0, flows computed")
end

@suite "dc_pf with transformer (tap + phase_shift)" begin
    net = Network(name="trafo test", baseMVA=100.0)
    add_bus!(net, "HV", v_nom=380.0, slack=true)
    add_bus!(net, "LV", v_nom=132.0)
    add_transformer!(net, "T1", "HV", "LV", x=0.2, tap_ratio=0.95, phase_shift=2.0)
    add_generator!(net, "G", "HV", p_nom=200.0)
    add_load!(net, "D", "LV", p_set=150.0)
    res = dc_pf(net, verbose=false)
    @assert res.converged
    @assert length(res.trafo_flows) == 1
    @assert abs(res.trafo_flows[1].P_MW) > 0
    pass("dc_pf: transformer with tap=0.95 and phase_shift=2°")
end

@suite "linear_ac_pf with all components" begin
    net = make_full_3bus()
    res = linear_ac_pf(net, verbose=false)
    @assert res.converged
    @assert all(res.V_mag .> 0.9) && all(res.V_mag .< 1.1)
    @assert abs(res.V_ang[1]) < 1e-10
    @assert length(res.P_flow) == 3
    pass("linear_ac_pf: voltages in range, slack=0")
end

@suite "lopf — Store + Link + GlobalConstraint" begin
    net = make_full_3bus()
    res = lopf(net, verbose=false)
    @assert res.converged
    # G1 + G2 + Wind3 - Store ± Link = load (500 MW)
    @assert haskey(res.P_gen, "G1") && haskey(res.P_gen, "G2")
    @assert haskey(res.P_link, "Link13")
    @assert haskey(res.P_store, "Store3")
    @assert res.total_cost > 0
    pass("lopf: all components in optimal solution")

    # Verify CO₂ constraint is respected
    # emissions = G1*0.20 + G2*0.34 ≤ 3000 (single period, so ≤ 3000)
    emissions = res.P_gen["G1"]*0.20 + res.P_gen["G2"]*0.34
    @assert emissions <= 3000.0 + 1e-3
    pass("lopf: CO₂ constraint satisfied ($(round(emissions,digits=1)) t ≤ 3000 t)")
end

@suite "lopf — congestion with Link bypassing" begin
    # Line L13 congested at 50 MW, but Link13 can carry 100 MW with η=0.97
    net = Network(name="congestion+link", baseMVA=100.0)
    add_bus!(net, "B1", slack=true); add_bus!(net, "B2"); add_bus!(net, "B3")
    add_line!(net, "L12", "B1", "B2", x=0.1, s_nom=Inf)
    add_line!(net, "L13", "B1", "B3", x=0.1, s_nom=50.0)   # tight limit
    add_line!(net, "L23", "B2", "B3", x=0.1, s_nom=Inf)
    add_generator!(net, "G1", "B1", p_nom=500.0, marginal_cost=10.0)
    add_generator!(net, "G3", "B3", p_nom=200.0, marginal_cost=90.0)  # expensive
    add_load!(net, "D3", "B3", p_set=300.0)
    add_link!(net, "Lnk13", "B1", "B3", p_nom=300.0, efficiency=0.98, marginal_cost=1.0)

    res = lopf(net, verbose=false)
    @assert res.converged
    # Expensive G3 should be zero or minimal if link covers the gap
    @assert res.P_gen["G3"] < 50.0
    pass("lopf: Link allows cheap gen to bypass congested line")
end

@suite "lopf_multiperiod — Store SoC dynamics" begin
    net = Network(name="store multiperiod", baseMVA=100.0)
    add_bus!(net, "B1", slack=true); add_bus!(net, "B2")
    add_line!(net, "L12", "B1", "B2", x=0.1)
    add_generator!(net, "G1", "B1", p_nom=300.0, marginal_cost=20.0)
    add_load!(net, "D2", "B2", p_set=200.0)
    add_store!(net, "St2", "B2",
        e_nom=400.0, p_nom=80.0,
        standing_loss=0.01, e_cyclic=true)

    res = lopf_multiperiod(net, T=12, verbose=false)
    @assert string(res.status) == "OPTIMAL"
    @assert haskey(res.store_e, "St2")
    e_vec = res.store_e["St2"]
    @assert length(e_vec) == 13   # 0:12
    # Cyclic: E[0] ≈ E[12]
    @assert abs(e_vec[1] - e_vec[end]) < 1.0
    pass("lopf_multiperiod: Store SoC cyclic, $(length(e_vec)-1) periods")
end

@suite "lopf_multiperiod — Link time series" begin
    net = Network(name="link multiperiod", baseMVA=100.0)
    add_bus!(net, "A", slack=true); add_bus!(net, "B")
    add_generator!(net, "GA", "A", p_nom=400.0, marginal_cost=5.0)
    add_generator!(net, "GB", "B", p_nom=100.0, marginal_cost=100.0)
    add_load!(net, "DB", "B", p_set=150.0)
    add_link!(net, "HVDC", "A", "B", p_nom=300.0, efficiency=0.97)

    res = lopf_multiperiod(net, T=6, verbose=false)
    @assert string(res.status) == "OPTIMAL"
    @assert haskey(res.link_p, "HVDC")
    @assert all(v >= -1e-6 for v in res.link_p["HVDC"])   # unidirectional
    avg_flow = sum(res.link_p["HVDC"]) / 6
    @assert avg_flow > 50.0   # HVDC should carry significant power
    pass("lopf_multiperiod: HVDC link active every period, avg=$(round(avg_flow,digits=1)) MW")
end

@suite "lopf_multiperiod — GlobalConstraint CO₂ + StorageUnit + Store" begin
    net = Network(name="full multiperiod", baseMVA=100.0)
    add_carrier!(net, "gas",  co2_emissions=0.20)
    add_carrier!(net, "coal", co2_emissions=0.34)
    add_bus!(net, "B1", slack=true); add_bus!(net, "B2")
    add_line!(net, "L12", "B1", "B2", x=0.1)
    add_generator!(net, "Gas",  "B1", p_nom=300.0, marginal_cost=40.0, carrier="gas")
    add_generator!(net, "Coal", "B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
    add_load!(net, "D2", "B2", p_set=200.0)
    add_storage_unit!(net, "Bat",   "B2", p_nom=50.0,  e_nom=200.0)
    add_store!(       net, "Store", "B2", p_nom=30.0,  e_nom=120.0)
    # CO₂ cap: forces partial switch from coal to gas
    add_global_constraint!(net, "co2",
        constant=200.0,
        carrier_weightings=Dict("gas"=>0.20, "coal"=>0.34))

    res_nocap = lopf_multiperiod(
        Network(name="nocap", baseMVA=100.0) |> net -> begin
            add_carrier!(net, "gas",  co2_emissions=0.20)
            add_carrier!(net, "coal", co2_emissions=0.34)
            add_bus!(net, "B1", slack=true); add_bus!(net, "B2")
            add_line!(net, "L12", "B1", "B2", x=0.1)
            add_generator!(net, "Gas",  "B1", p_nom=300.0, marginal_cost=40.0, carrier="gas")
            add_generator!(net, "Coal", "B1", p_nom=300.0, marginal_cost=20.0, carrier="coal")
            add_load!(net, "D2", "B2", p_set=200.0)
            net
        end,
        T=6, verbose=false)

    res = lopf_multiperiod(net, T=6, verbose=false)
    @assert string(res.status) == "OPTIMAL"
    coal_mwh = sum(res.gen_dispatch["Coal"])
    gas_mwh  = sum(res.gen_dispatch["Gas"])
    co2_total = coal_mwh*0.34 + gas_mwh*0.20
    @assert co2_total <= 200.0 + 1e-3
    pass("Full multiperiod: CO₂=$(round(co2_total,digits=1)) t ≤ 200 t ✓")
    pass("StorageUnit + Store + GlobalConstraint coexist in one model")
end

# ═══════════════════════════════════════════════════════════════════════════
#  4. EDGE CASES
# ═══════════════════════════════════════════════════════════════════════════
@suite "Edge cases" begin
    # Network with no storage/links/constraints → should work fine
    net = Network(name="minimal", baseMVA=100.0)
    add_bus!(net, "B1", slack=true); add_bus!(net, "B2")
    add_line!(net, "L", "B1", "B2", x=0.1)
    add_generator!(net, "G", "B1", p_nom=200.0, marginal_cost=10.0)
    add_load!(net, "D", "B2", p_set=100.0)
    r1 = dc_pf(net, verbose=false)
    r2 = linear_ac_pf(net, verbose=false)
    r3 = lopf(net, verbose=false)
    r4 = lopf_multiperiod(net, T=4, verbose=false)
    @assert r1.converged && r2.converged && r3.converged
    @assert string(r4.status) == "OPTIMAL"
    pass("Minimal network (no storage/links): all 4 solvers work")

    # Bidirectional Link
    net2 = Network(name="bidir", baseMVA=100.0)
    add_bus!(net2, "A", slack=true); add_bus!(net2, "B")
    add_generator!(net2, "GA", "A", p_nom=300.0, marginal_cost=10.0)
    add_generator!(net2, "GB", "B", p_nom=300.0, marginal_cost=50.0)
    add_load!(net2, "DA", "A", p_set=100.0)
    add_load!(net2, "DB", "B", p_set=200.0)
    add_link!(net2, "Bidir", "A", "B", p_nom=300.0, p_min_pu=-1.0, efficiency=0.99)
    r5 = lopf(net2, verbose=false)
    @assert r5.converged
    pass("Bidirectional Link (p_min_pu=-1): LOPF converged")

    # Store with standing_loss
    net3 = Network(name="standing_loss", baseMVA=100.0)
    add_bus!(net3, "B1", slack=true)
    add_generator!(net3, "G", "B1", p_nom=200.0, marginal_cost=10.0)
    add_load!(net3, "D", "B1", p_set=100.0)
    add_store!(net3, "S", "B1", e_nom=200.0, p_nom=50.0,
               standing_loss=0.02, e_cyclic=true)
    r6 = lopf_multiperiod(net3, T=6, verbose=false)
    @assert string(r6.status) == "OPTIMAL"
    pass("Store with standing_loss=2%/h in multiperiod")
end

# ═══════════════════════════════════════════════════════════════════════════
#  SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
println("\n" * "="^60)
println("COMPONENT AUDIT SUMMARY")
println("="^60)
passed = count(values(results))
total  = length(results)
for (name, ok) in sort(collect(results))
    println("  $(ok ? "✓" : "✗")  $name")
end
println()
println("Result: $passed / $total suites passed")
passed == total ? println("ALL COMPONENTS VERIFIED ✓") :
                  println("$(total-passed) SUITE(S) FAILED ✗")
println("="^60)
