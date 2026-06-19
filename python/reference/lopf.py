"""
Linear Optimal Power Flow (LOPF) via PyPSA
Для сравнения с Julia-реализацией (JuMP + HiGHS)
"""
import pypsa
import numpy as np

def solve_lopf(line_capacity=None, verbose=True):
    """
    3-узловая сеть с 2 генераторами разной стоимости.

    Сценарий А: без ограничений линий → всё от G1 (дешёвый)
    Сценарий Б: ограничение 150 МВт → G2 вынужден работать
    """
    n = pypsa.Network()

    # Узлы
    n.add("Bus", "Bus 1", v_nom=380.0)
    n.add("Bus", "Bus 2", v_nom=380.0)
    n.add("Bus", "Bus 3", v_nom=380.0)

    # Линии
    s_nom = line_capacity if line_capacity is not None else 1e6   # МВА
    n.add("Line", "Line 1-2", bus0="Bus 1", bus1="Bus 2",
          x=0.1, r=0.01, s_nom=s_nom)
    n.add("Line", "Line 1-3", bus0="Bus 1", bus1="Bus 3",
          x=0.1, r=0.01, s_nom=s_nom)
    n.add("Line", "Line 2-3", bus0="Bus 2", bus1="Bus 3",
          x=0.1, r=0.01, s_nom=s_nom)

    # Генераторы (с marginal_cost)
    n.add("Generator", "G1", bus="Bus 1",
          p_nom=400.0, marginal_cost=20.0, control="Slack")
    n.add("Generator", "G2", bus="Bus 2",
          p_nom=300.0, marginal_cost=50.0)

    # Нагрузки (снэпшот = 1 час)
    n.add("Load", "Load 2", bus="Bus 2", p_set=200.0)
    n.add("Load", "Load 3", bus="Bus 3", p_set=300.0)

    # Решение LOPF
    n.optimize(solver_name="highs")

    if verbose:
        print("=" * 60)
        print("LOPF RESULTS (PyPSA)")
        print("=" * 60)

        print("\n1. GENERATOR DISPATCH:")
        print(f"{'Generator':<20} {'P (MW)':>10} {'P_max (MW)':>12} {'Cost (€/MWh)':>14}")
        print("-" * 58)
        for gen_name in n.generators.index:
            p    = n.generators_t.p[gen_name].values[0]
            pmax = n.generators.loc[gen_name, "p_nom"]
            cost = n.generators.loc[gen_name, "marginal_cost"]
            print(f"{gen_name:<20} {p:>10.2f} {pmax:>12.2f} {cost:>14.2f}")

        print("\n2. LINE FLOWS:")
        print(f"{'Line':<15} {'P (MW)':>10} {'P_max (MW)':>12} {'Loading (%)':>12}")
        print("-" * 52)
        for line_name in n.lines.index:
            p_flow = n.lines_t.p0[line_name].values[0]
            s_nom_line = n.lines.loc[line_name, "s_nom"]
            loading = abs(p_flow) / s_nom_line * 100 if s_nom_line < 1e5 else float("nan")
            print(f"{line_name:<15} {p_flow:>10.2f} "
                  f"{'∞' if s_nom_line > 1e5 else f'{s_nom_line:.1f}':>12} "
                  f"{'—' if np.isnan(loading) else f'{loading:.1f}%':>12}")

        total_cost = (n.generators_t.p * n.generators.marginal_cost).sum(axis=1).sum()
        print(f"\nTotal generation cost: {total_cost:.2f} €/h")
        print("=" * 60)

    # Возвращаем результаты для сравнения
    p_gen = {gen: n.generators_t.p[gen].values[0]
             for gen in n.generators.index}
    p_line = {line: n.lines_t.p0[line].values[0]
              for line in n.lines.index}
    total_cost = (n.generators_t.p * n.generators.marginal_cost).sum(axis=1).sum()

    return {"p_gen": p_gen, "p_line": p_line, "total_cost": total_cost}


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("LINEAR OPTIMAL POWER FLOW (LOPF)")
    print("Python/PyPSA implementation")
    print("=" * 60)

    # Сценарий А: без ограничений
    print("\n[A] SCENARIO A: No line capacity limits")
    print("   (Expected: G1 = 400 MW, G2 = 100 MW, cost = 13,000 e/h)")
    res_A = solve_lopf(line_capacity=None)

    # Сценарий Б: ограничение 200 МВт
    print("\n[B] SCENARIO B: Line capacity = 200 MW")
    print("   (G1 constrained by grid -> G2 must compensate)")
    res_B = solve_lopf(line_capacity=200.0)

    # Сравнение сценариев
    print("\n" + "=" * 60)
    print("COMPARISON: Unconstrained vs Constrained")
    print("=" * 60)
    print(f"\n{'Metric':<25} {'Scenario A':>15} {'Scenario B':>15}")
    print("-" * 57)
    print(f"{'G1 dispatch (MW)':<25} {res_A['p_gen']['G1']:>15.2f} {res_B['p_gen']['G1']:>15.2f}")
    print(f"{'G2 dispatch (MW)':<25} {res_A['p_gen']['G2']:>15.2f} {res_B['p_gen']['G2']:>15.2f}")
    print(f"{'Total cost (€/h)':<25} {res_A['total_cost']:>15.2f} {res_B['total_cost']:>15.2f}")
    increase = (res_B['total_cost'] - res_A['total_cost']) / res_A['total_cost'] * 100
    print(f"{'Cost increase':<25} {'—':>15} {increase:>14.1f}%")
    print(f"\n[OK] Congestion increases cost by {increase:.1f}%!")
