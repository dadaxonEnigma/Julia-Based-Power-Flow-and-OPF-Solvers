"""
Linearized AC Power Flow — Python reference implementation.
Validates against PyPSA full AC PF (network.pf()).
"""
import numpy as np
import pypsa
import warnings
import logging

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)


def solve_linear_ac_pf(buses, lines, generators, loads, baseMVA=100.0, v_nom=380.0):
    """
    Linearized AC Power Flow around flat start (|V|=1, theta=0).

    Solves:
        [ B' -G'] [Δθ ]   [P]
        [-G' -B'] [Δ|V|] = [Q]

    Parameters
    ----------
    buses      : list of bus names
    lines      : list of (from, to, r_pu, x_pu)
    generators : dict {bus: (P_MW, Q_MVAr)}
    loads      : dict {bus: (P_MW, Q_MVAr)}
    baseMVA    : system base [MVA]
    """
    n = len(buses)
    # r, x in physical Ω → convert to p.u.
    z_base = v_nom**2 / baseMVA    # 380²/100 = 1444 Ω
    G = np.zeros((n, n))
    B = np.zeros((n, n))

    for (f, t, r, x) in lines:
        f -= 1; t -= 1
        r_pu  = r / z_base
        x_pu  = x / z_base
        denom = r_pu**2 + x_pu**2
        g_ij  =  r_pu / denom
        b_ij  =  x_pu / denom      # positive |susceptance|

        G[f,f] += g_ij;  G[t,t] += g_ij
        B[f,f] += b_ij;  B[t,t] += b_ij
        G[f,t] -= g_ij;  G[t,f] -= g_ij
        B[f,t] -= b_ij;  B[t,f] -= b_ij

    G_mw = G * baseMVA
    B_mw = B * baseMVA

    P = np.zeros(n)
    Q = np.zeros(n)
    for bus, (p, q) in generators.items():
        P[bus-1] += p;  Q[bus-1] += q
    for bus, (p, q) in loads.items():
        P[bus-1] -= p;  Q[bus-1] -= q

    # Remove slack bus (index 0)
    idx = list(range(1, n))
    B_r = B_mw[np.ix_(idx, idx)]
    G_r = G_mw[np.ix_(idx, idx)]
    P_r = P[idx]
    Q_r = Q[idx]

    nm = n - 1
    A = np.block([[ B_r, -G_r],
                  [-G_r, -B_r]])
    rhs = np.concatenate([P_r, Q_r])
    x_sol = np.linalg.solve(A, rhs)

    dtheta = np.zeros(n);  dtheta[idx] = x_sol[:nm]
    dv     = np.zeros(n);  dv[idx]     = x_sol[nm:]
    V_mag  = 1.0 + dv
    theta  = dtheta

    P_flow, Q_flow = [], []
    for (f, t, r, x) in lines:
        f -= 1; t -= 1
        r_pu  = r / z_base
        x_pu  = x / z_base
        denom  = r_pu**2 + x_pu**2
        g_ij   = r_pu / denom
        b_ij   = x_pu / denom      # positive |susceptance|
        dth_ft = theta[f] - theta[t]
        dv_ft  = dv[f] - dv[t]
        P_flow.append(( b_ij * dth_ft + g_ij * dv_ft) * baseMVA)
        Q_flow.append(( b_ij * dv_ft  - g_ij * dth_ft) * baseMVA)

    return V_mag, theta, np.array(P_flow), np.array(Q_flow)


# ----------------------------------------------------------------
#  Full AC PF reference via PyPSA
# ----------------------------------------------------------------
def pypsa_ac_pf():
    net = pypsa.Network()
    net.set_snapshots([0])
    for i, name in enumerate(["Bus 0","Bus 1","Bus 2"]):
        net.add("Bus", name, v_nom=380.0)
    for k,(f,t,r,x) in enumerate([(1,2,0.01,0.1),(1,3,0.01,0.1),(2,3,0.01,0.1)]):
        net.add("Line", f"L{k}", bus0=f"Bus {f-1}", bus1=f"Bus {t-1}",
                r=r, x=x, s_nom=1e6)
    net.add("Generator","G0",bus="Bus 0",p_nom=500,control="Slack")
    net.add("Load","L1",bus="Bus 1",p_set=300)
    net.add("Load","L2",bus="Bus 2",p_set=200)
    net.pf()
    buses = ["Bus 0","Bus 1","Bus 2"]
    vm  = net.buses_t.v_mag_pu.loc[0, buses].values
    va  = net.buses_t.v_ang.loc[0, buses].values
    pf  = net.lines_t.p0.loc[0].values
    qf  = net.lines_t.q0.loc[0].values
    return vm, va, pf, qf


# ================================================================
SEP = "=" * 60
print(SEP)
print("LINEARIZED AC POWER FLOW  (Python reference)")
print(SEP)

buses = ["Bus 0","Bus 1","Bus 2"]
lines = [(1,2,0.01,0.1),(1,3,0.01,0.1),(2,3,0.01,0.1)]
generators = {1: (500.0, 0.0)}
loads      = {2: (300.0, 0.0), 3: (200.0, 0.0)}

vm_lac, va_lac, pf_lac, qf_lac = solve_linear_ac_pf(
    buses, lines, generators, loads)

print("\nRunning PyPSA full AC PF for reference...")
vm_ref, va_ref, pf_ref, qf_ref = pypsa_ac_pf()

print(f"\n{'Bus':<8} {'|V| LACPF':>10} {'|V| Full AC':>11} {'|Δ|':>10}")
print("-"*42)
for i in range(3):
    print(f"{buses[i]:<8} {vm_lac[i]:>10.6f} {vm_ref[i]:>11.6f} "
          f"{abs(vm_lac[i]-vm_ref[i]):>10.2e}")

print(f"\n{'Bus':<8} {'θ LACPF (rad)':>13} {'θ Full AC (rad)':>15} {'|Δ|':>10}")
print("-"*50)
for i in range(3):
    print(f"{buses[i]:<8} {va_lac[i]:>13.6f} {va_ref[i]:>15.6f} "
          f"{abs(va_lac[i]-va_ref[i]):>10.2e}")

lnames = ["L 0-1","L 0-2","L 1-2"]
print(f"\n{'Line':<8} {'P LACPF':>10} {'P Full AC':>10} {'|ΔP|':>8} "
      f"{'Q LACPF':>10} {'Q Full AC':>10} {'|ΔQ|':>8}")
print("-"*68)
for i in range(3):
    print(f"{lnames[i]:<8} {pf_lac[i]:>10.4f} {pf_ref[i]:>10.4f} "
          f"{abs(pf_lac[i]-pf_ref[i]):>8.4f} "
          f"{qf_lac[i]:>10.4f} {qf_ref[i]:>10.4f} "
          f"{abs(qf_lac[i]-qf_ref[i]):>8.4f}")

print(f"\nMax |ΔVm| = {max(abs(vm_lac-vm_ref)):.2e} p.u.")
print(f"Max |Δθ|  = {max(abs(va_lac-va_ref)):.2e} rad")
print(f"Max |ΔP|  = {max(abs(pf_lac-pf_ref)):.4f} MW")
print(f"Max |ΔQ|  = {max(abs(qf_lac-qf_ref)):.4f} MVAr")
print(SEP)
