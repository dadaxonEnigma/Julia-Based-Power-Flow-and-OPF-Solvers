using Plots
using GraphPlot
using LinearAlgebra
using Printf

include("../solvers/dc_power_flow.jl")

"""
Визуализация результатов DC Power Flow
"""
function visualize_network_results(network::PowerNetwork, results)
    
    # 1. ТОПОЛОГИЯ СЕТИ
    println("\n" * "="^60)
    println("NETWORK TOPOLOGY VISUALIZATION")
    println("="^60)
    
    # Создаем граф из линий
    n = network.n_buses
    adj_matrix = zeros(Int, n, n)
    
    for (from, to, x) in network.lines
        adj_matrix[from, to] = 1
        adj_matrix[to, from] = 1
    end
    
    # Позиции узлов (простая круговая раскладка)
    θ_layout = range(0, 2π, length=n+1)[1:n]
    x_pos = cos.(θ_layout)
    y_pos = sin.(θ_layout)
    
    # График 1: Топология сети
    p1 = plot(size=(800, 800), title="Network Topology", 
              legend=:topright, aspect_ratio=:equal)
    
    # Рисуем линии
    for (idx, (from, to, x)) in enumerate(network.lines)
        plot!(p1, [x_pos[from], x_pos[to]], [y_pos[from], y_pos[to]], 
              color=:gray, linewidth=2, label="")
    end
    
    # Рисуем узлы
    scatter!(p1, x_pos, y_pos, 
            markersize=20, 
            markercolor=:lightblue,
            markerstrokewidth=2,
            markerstrokecolor=:black,
            label="Buses")
    
    # Подписи узлов
    for i in 1:n
        annotate!(p1, x_pos[i], y_pos[i], 
                 text(network.buses[i], 8, :center))
    end
    
    # Отмечаем генераторы
    for (bus, power) in network.generators
        scatter!(p1, [x_pos[bus]], [y_pos[bus]], 
                markersize=25, 
                markercolor=:green,
                markershape=:star5,
                label=bus==first(keys(network.generators)) ? "Generators" : "")
    end
    
    # Отмечаем нагрузки
    for (bus, power) in network.loads
        scatter!(p1, [x_pos[bus]], [y_pos[bus]], 
                markersize=25, 
                markercolor=:red,
                markershape=:diamond,
                label=bus==first(keys(network.loads)) ? "Loads" : "")
    end
    
    display(p1)
    savefig(p1, "network_topology.png")
    println("✓ Saved: network_topology.png")
    
    # 2. УГЛЫ НАПРЯЖЕНИЯ
    println("\n" * "="^60)
    println("VOLTAGE ANGLES VISUALIZATION")
    println("="^60)
    
    p2 = bar(1:n, results.θ .* (180/π), 
            xlabel="Bus", 
            ylabel="Voltage Angle (degrees)",
            title="Voltage Angles at Each Bus",
            legend=false,
            color=:blue,
            xticks=(1:n, network.buses),
            xrotation=45)
    
    hline!(p2, [0], color=:red, linestyle=:dash, linewidth=2)
    
    display(p2)
    savefig(p2, "voltage_angles.png")
    println("✓ Saved: voltage_angles.png")
    
    # 3. ПОТОКИ НА ЛИНИЯХ
    println("\n" * "="^60)
    println("LINE FLOWS VISUALIZATION")
    println("="^60)
    
    p3 = bar(1:length(network.line_names), results.line_flows,
            xlabel="Line",
            ylabel="Power Flow (MW)",
            title="Power Flow on Transmission Lines",
            legend=false,
            color=ifelse.(results.line_flows .>= 0, :green, :red),
            xticks=(1:length(network.line_names), network.line_names),
            xrotation=45)
    
    hline!(p3, [0], color=:black, linewidth=1)
    
    display(p3)
    savefig(p3, "line_flows.png")
    println("✓ Saved: line_flows.png")
    
    # 4. СЕТЬ С ПОТОКАМИ
    println("\n" * "="^60)
    println("NETWORK WITH FLOWS")
    println("="^60)
    
    p4 = plot(size=(800, 800), 
              title="Network with Power Flows", 
              legend=:topright, 
              aspect_ratio=:equal)
    
    # Нормализуем потоки для визуализации
    max_flow = maximum(abs.(results.line_flows))
    
    # Рисуем линии с толщиной пропорциональной потоку
    for (idx, (from, to, x)) in enumerate(network.lines)
        flow = results.line_flows[idx]
        linewidth = 1 + 5 * abs(flow) / max_flow
        color = flow >= 0 ? :green : :red
        
        plot!(p4, [x_pos[from], x_pos[to]], [y_pos[from], y_pos[to]], 
              color=color, linewidth=linewidth, 
              label="", arrow=true)
        
        # Подпись потока на линии
        mid_x = (x_pos[from] + x_pos[to]) / 2
        mid_y = (y_pos[from] + y_pos[to]) / 2
        annotate!(p4, mid_x, mid_y, 
                 text(@sprintf("%.1f", flow), 7, :center))
    end
    
    # Рисуем узлы с цветом по углу
    angles_norm = (results.θ .- minimum(results.θ)) ./ (maximum(results.θ) - minimum(results.θ) .+ 1e-10)
    
    scatter!(p4, x_pos, y_pos, 
            markersize=20, 
            marker_z=results.θ .* (180/π),
            color=:viridis,
            markerstrokewidth=2,
            markerstrokecolor=:black,
            label="Bus angles",
            colorbar_title="Angle (°)")
    
    # Подписи узлов
    for i in 1:n
        annotate!(p4, x_pos[i], y_pos[i] + 0.15, 
                 text(network.buses[i], 8, :center))
    end
    
    display(p4)
    savefig(p4, "network_with_flows.png")
    println("✓ Saved: network_with_flows.png")
    
    # 5. СВОДНАЯ ТАБЛИЦА
    println("\n" * "="^60)
    println("SUMMARY TABLE")
    println("="^60)
    
    println("\nBus Summary:")
    println("="^50)
    @printf("%-10s %12s %12s %12s\n", "Bus", "Angle(°)", "Gen(MW)", "Load(MW)")
    println("="^50)
    
    for i in 1:n
        gen = get(network.generators, i, 0.0)
        load = get(network.loads, i, 0.0)
        @printf("%-10s %12.6f %12.1f %12.1f\n", 
                network.buses[i], 
                results.θ[i] * 180/π,
                gen,
                load)
    end
    
    println("\nLine Summary:")
    println("="^50)
    @printf("%-15s %12s %12s\n", "Line", "Flow(MW)", "Direction")
    println("="^50)
    
    for (idx, name) in enumerate(network.line_names)
        flow = results.line_flows[idx]
        direction = flow >= 0 ? "→" : "←"
        @printf("%-15s %12.2f %12s\n", name, abs(flow), direction)
    end
    
    println("\n" * "="^60)
    println("ALL VISUALIZATIONS COMPLETED!")
    println("="^60)
end

# Запуск с нашей тестовой сетью
println("\n\nRunning DC Power Flow with Visualization...")

# Та же сеть что и раньше
buses = ["Bus 0", "Bus 1", "Bus 2"]
lines = [
    (1, 2, 0.1),
    (1, 3, 0.1),
    (2, 3, 0.1)
]
line_names = ["Line 0-1", "Line 0-2", "Line 1-2"]
generators = Dict(1 => 500.0)
loads = Dict(2 => 300.0, 3 => 200.0)

network = PowerNetwork(3, buses, lines, line_names, generators, loads, 1)

# Решаем
results = solve_dc_power_flow(network)

# Визуализируем
visualize_network_results(network, results)