"""
visualize_results.jl  —  Thesis figures for DC PF + LOPF validation results
Generates ../results/results_validation.png  (1×3 panel, 300 DPI, improved labels)
"""

using Plots
using StatsPlots
using Printf

# ─── 3-bus test network ───────────────────────────────────────────────────────
b    = 1000.0
B    = [2b  -b  -b;  -b  2b  -b;  -b  -b  2b]
Bred = B[2:3, 2:3]

# ── DC Power Flow ─────────────────────────────────────────────────────────────
P_dcpf = [400.0, 100.0 - 200.0, 0.0 - 300.0]
θ_dc   = vcat(0.0, Bred \ P_dcpf[2:3])
pf_dc  = [b*(θ_dc[1]-θ_dc[2]), b*(θ_dc[1]-θ_dc[3]), b*(θ_dc[2]-θ_dc[3])]

# ── LOPF scenarios ────────────────────────────────────────────────────────────
G1_A, G2_A, cost_A = 400.0, 100.0, 13_000.0
pf_A = pf_dc
P_B  = [300.0, 200.0 - 200.0, 0.0 - 300.0]
θ_B  = vcat(0.0, Bred \ P_B[2:3])
pf_B = [b*(θ_B[1]-θ_B[2]), b*(θ_B[1]-θ_B[3]), b*(θ_B[2]-θ_B[3])]
G1_B, G2_B, cost_B = 300.0, 200.0, 16_000.0
pct_increase = (cost_B - cost_A) / cost_A * 100

# ─── Colors ───────────────────────────────────────────────────────────────────
C_GEN  = RGB(0.13, 0.55, 0.13)
C_BOTH = RGB(0.85, 0.50, 0.10)
C_LOAD = RGB(0.80, 0.15, 0.15)
C_A    = RGB(0.13, 0.47, 0.71)
C_B    = RGB(0.80, 0.15, 0.15)
DPI    = 300

BASE = (tickfontsize=11, guidefontsize=13, titlefontsize=14,
        legendfontsize=11, framestyle=:box, grid=true, gridalpha=0.25,
        foreground_color_legend=nothing,
        background_color_legend=RGBA(1,1,1,0.85))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PANEL A — Voltage angles (DC Power Flow)                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
θ_deg = θ_dc .* (180 / π)

pA = bar(1:3, θ_deg;
    color  = [C_GEN, C_BOTH, C_LOAD],
    alpha  = 0.82,
    label  = "",
    title  = "DC Power Flow\nVoltage Angles",
    xlabel = "Bus",
    ylabel = "Angle  θ  (degrees)",
    xticks = (1:3, ["Bus 1", "Bus 2", "Bus 3"]),
    ylims  = (-17, 4),
    BASE..., dpi=DPI, size=(420, 480))

hline!(pA, [0.0]; color=:gray40, linestyle=:dash, linewidth=1.5,
       label="Reference  (0°)")

# Value annotations
for (i, θ) in enumerate(θ_deg)
    yoff = θ < 0 ? -1.5 : 1.0
    annotate!(pA, i, θ + yoff,
              Plots.text(@sprintf("%.2f°", θ), 11, :center))
end

# Bus type legend below
annotate!(pA, 1, -16, Plots.text("(slack, G₁)", 9, :center, C_GEN))
annotate!(pA, 2, -16, Plots.text("(G₂ + Load)", 9, :center, C_BOTH))
annotate!(pA, 3, -16, Plots.text("(Load)", 9, :center, C_LOAD))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PANEL B — LOPF generator dispatch: Scenario A vs B                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
dispatch_data = [G1_A  G1_B;
                 G2_A  G2_B]

pB = groupedbar(dispatch_data;
    bar_position = :dodge,
    bar_width    = 0.6,
    label        = ["Scenario A  (no limits)" "Scenario B  (limit 200 MW)"],
    color        = [C_A  C_B],
    alpha        = 0.82,
    title        = "LOPF\nGenerator Dispatch",
    xlabel       = "Generator",
    ylabel       = "Active power  (MW)",
    xticks       = (1:2, ["G₁", "G₂"]),
    ylims        = (0, 520),
    BASE..., dpi=DPI, size=(420, 480))

hline!(pB, [400.0]; linestyle=:dot, linewidth=1.5, color=C_GEN, label="")
hline!(pB, [300.0]; linestyle=:dot, linewidth=1.5, color=C_BOTH, label="")

# Value annotations above bars
for (g, vals) in enumerate([[G1_A, G1_B], [G2_A, G2_B]])
    for (s, val) in enumerate(vals)
        xoff = s == 1 ? -0.18 : 0.18
        yoff = 16
        annotate!(pB, g + xoff, val + yoff,
                  Plots.text(@sprintf("%.0f", val), 11, :center, :black))
    end
end

# Generator info legend
annotate!(pB, 1.5, 490, Plots.text("G₁: 20 €/MWh  |  G₂: 50 €/MWh", 9, :center))

# Cost annotations
annotate!(pB, 2.7, 450, Plots.text("Cost A:", 9, :right, C_A))
annotate!(pB, 2.7, 420, Plots.text("13 000 €/h", 9, :right, C_A))
annotate!(pB, 2.7, 380, Plots.text("Cost B:", 9, :right, C_B))
annotate!(pB, 2.7, 350, Plots.text(@sprintf("16 000 €/h  (+%.1f%%)", pct_increase), 9, :right, C_B))

# Capacity lines legend
annotate!(pB, 0.5, 405, Plots.text("capacity", 8, :center, :gray50))
annotate!(pB, 0.5, 295, Plots.text("capacity", 8, :center, :gray50))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PANEL C — LOPF line flows: Scenario A vs B                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
flow_data = [abs(pf_A[1])  abs(pf_B[1]);
             abs(pf_A[2])  abs(pf_B[2]);
             abs(pf_A[3])  abs(pf_B[3])]

pC = groupedbar(flow_data;
    bar_position = :dodge,
    bar_width    = 0.6,
    label        = ["Scenario A  (no limits)" "Scenario B  (limit 200 MW)"],
    color        = [C_A  C_B],
    alpha        = 0.82,
    title        = "LOPF\nLine Power Flows",
    xlabel       = "Line",
    ylabel       = "|Power flow|  (MW)",
    xticks       = (1:3, ["1-2", "1-3", "2-3"]),
    ylims        = (0, 320),
    BASE..., dpi=DPI, size=(420, 480))

hline!(pC, [200.0]; linestyle=:dash, linewidth=2, color=:gray30,
       label="Thermal limit")

# Value annotations above bars
for (ln, (fA, fB)) in enumerate(zip(abs.(pf_A), abs.(pf_B)))
    annotate!(pC, ln - 0.18, fA + 10,
              Plots.text(@sprintf("%.0f", fA), 11, :center, C_A))
    annotate!(pC, ln + 0.18, fB + 10,
              Plots.text(@sprintf("%.0f", fB), 11, :center, C_B))
end

# Annotations for line labels
annotate!(pC, 1, -28, Plots.text("Line 1-2", 9, :center))
annotate!(pC, 2, -28, Plots.text("Line 1-3", 9, :center))
annotate!(pC, 3, -28, Plots.text("Line 2-3", 9, :center))

# Mark congested line
annotate!(pC, 2.7, 245, Plots.text("⚠ Congested", 10, C_B, :right))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Combined 1×3 figure                                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
pfig = plot(pA, pB, pC;
    layout        = (1, 3),
    size          = (1450, 540),
    plot_title    = "Validation: DC Power Flow & LOPF on 3-Bus Test Network",
    plot_titlefontsize = 15,
    left_margin   = 10Plots.mm,
    right_margin  = 8Plots.mm,
    bottom_margin = 12Plots.mm,
    top_margin    = 8Plots.mm,
    dpi           = DPI)

savefig(pfig, "../results/results_validation.png")
println("Saved: results/results_validation.png  (improved labels)")

# ─── Summary ──────────────────────────────────────────────────────────────────
println("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
println("  Validation Results Summary")
println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
println("\n  DC Power Flow  (G1=400 MW, G2=100 MW):")
@printf("    θ₁ = 0.00°  (slack)\n")
@printf("    θ₂ = %.2f°\n", θ_dc[2]*180/π)
@printf("    θ₃ = %.2f°\n", θ_dc[3]*180/π)
@printf("    p₁₂ = %.1f MW  |  p₁₃ = %.1f MW  |  p₂₃ = %.1f MW\n", pf_dc...)
println("\n  LOPF Scenario A  (Unconstrained):")
@printf("    G1 = %.0f MW  |  G2 = %.0f MW\n", G1_A, G2_A)
@printf("    Cost = %.0f €/h\n", cost_A)
@printf("    Line 1-3: %.1f MW  (%.1f%% of 200 MW limit)\n", abs(pf_A[2]), abs(pf_A[2])/2)
println("\n  LOPF Scenario B  (200 MW Thermal Limit):")
@printf("    G1 = %.0f MW  |  G2 = %.0f MW\n", G1_B, G2_B)
@printf("    Cost = %.0f €/h  (+%.1f%%)\n", cost_B, pct_increase)
@printf("    Line 1-3: %.1f MW  (%.1f%% of 200 MW limit) ← at limit\n", abs(pf_B[2]), abs(pf_B[2])/2)
println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
