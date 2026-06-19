import pypsa
import numpy as np
import pandas as pd

print("="*60)
print("TESTING PyPSA - DC POWER FLOW")
print("="*60)

# Создаем простую сеть: 3 узла, 3 линии
network = pypsa.Network()

# Добавляем узлы (buses)
network.add("Bus", "Bus 0", v_nom=380)
network.add("Bus", "Bus 1", v_nom=380)  
network.add("Bus", "Bus 2", v_nom=380)

print("\n✓ Added 3 buses")

# Добавляем линии
network.add("Line", "Line 0-1",
            bus0="Bus 0", bus1="Bus 1",
            x=0.1, r=0.01, s_nom=1000)

network.add("Line", "Line 0-2",
            bus0="Bus 0", bus1="Bus 2",
            x=0.1, r=0.01, s_nom=1000)

network.add("Line", "Line 1-2",
            bus0="Bus 1", bus1="Bus 2",
            x=0.1, r=0.01, s_nom=1000)

print("✓ Added 3 lines")

# Добавляем генератор на Bus 0
network.add("Generator", "Gen 0",
            bus="Bus 0",
            p_nom=500,
            marginal_cost=10)

print("✓ Added 1 generator")

# Добавляем нагрузки
network.add("Load", "Load 1",
            bus="Bus 1",
            p_set=300)

network.add("Load", "Load 2",
            bus="Bus 2",
            p_set=200)

print("✓ Added 2 loads")

print("\n" + "="*60)
print("NETWORK STRUCTURE")
print("="*60)
print(f"Buses: {len(network.buses)}")
print(f"Lines: {len(network.lines)}")
print(f"Generators: {len(network.generators)}")
print(f"Loads: {len(network.loads)}")

print("\n" + "="*60)
print("RUNNING DC POWER FLOW...")
print("="*60)

# Запускаем DC Power Flow
network.lpf()

print("\n✓ DC Power Flow completed successfully!")

print("\n" + "="*60)
print("RESULTS")
print("="*60)

print("\n1. BUS VOLTAGE ANGLES (radians):")
print(network.buses_t.v_ang)

print("\n2. LINE POWER FLOWS (MW):")
print(network.lines_t.p0)

print("\n3. GENERATOR OUTPUT (MW):")
print(network.generators_t.p)

print("\n4. SUMMARY:")
total_generation = network.generators_t.p.sum().sum()
total_load = network.loads.p_set.sum()
print(f"Total Generation: {total_generation:.2f} MW")
print(f"Total Load: {total_load:.2f} MW")
print(f"Balance: {total_generation - total_load:.2f} MW (should be ~0)")

print("\n" + "="*60)
print("TEST COMPLETED SUCCESSFULLY!")
print("="*60)