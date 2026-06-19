import pypsa
import numpy as np

print("="*60)
print("TESTING PyPSA - AC POWER FLOW")
print("="*60)

# Создаем ту же сеть
network = pypsa.Network()

# Узлы
network.add("Bus", "Bus 0", v_nom=380)
network.add("Bus", "Bus 1", v_nom=380)  
network.add("Bus", "Bus 2", v_nom=380)

# Линии (добавляем сопротивление r для AC PF)
network.add("Line", "Line 0-1",
            bus0="Bus 0", bus1="Bus 1",
            x=0.1, r=0.01, s_nom=1000)

network.add("Line", "Line 0-2",
            bus0="Bus 0", bus1="Bus 2",
            x=0.1, r=0.01, s_nom=1000)

network.add("Line", "Line 1-2",
            bus0="Bus 1", bus1="Bus 2",
            x=0.1, r=0.01, s_nom=1000)

# Генератор - указываем control type
network.add("Generator", "Gen 0",
            bus="Bus 0",
            p_nom=500,
            control="Slack")  # Slack bus для AC PF

# Нагрузки
network.add("Load", "Load 1",
            bus="Bus 1",
            p_set=300)

network.add("Load", "Load 2",
            bus="Bus 2",
            p_set=200)

print("\n" + "="*60)
print("RUNNING AC POWER FLOW...")
print("="*60)

# Запускаем AC Power Flow (не linear!)
network.pf()  # Это AC PF, не DC!

print("\n✓ AC Power Flow completed!")

print("\n" + "="*60)
print("RESULTS")
print("="*60)

print("\n1. BUS VOLTAGE MAGNITUDES (p.u.):")
print(network.buses_t.v_mag_pu)

print("\n2. BUS VOLTAGE ANGLES (radians):")
print(network.buses_t.v_ang)

print("\n3. LINE ACTIVE POWER FLOWS (MW):")
print(network.lines_t.p0)

print("\n4. LINE REACTIVE POWER FLOWS (MVAr):")
print(network.lines_t.q0)

print("\n5. GENERATOR OUTPUT:")
print("Active power (P):")
print(network.generators_t.p)
print("\nReactive power (Q):")
print(network.generators_t.q)

print("\n6. CONVERGENCE INFO:")
if hasattr(network, 'pf_info'):
    print(f"Iterations: {network.pf_info.get('n_iter', 'N/A')}")
    print(f"Converged: {network.pf_info.get('converged', 'N/A')}")

print("\n" + "="*60)
print("TEST COMPLETED!")
print("="*60)