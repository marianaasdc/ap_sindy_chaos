"""
hh_model_v2.py
==============
Hodgkin-Huxley model parameters, gate functions, ODE, and simulation
utilities shared by all HH-SINDy scripts.

Exports
-------
Parameters : gK, EK, gNa, ENa, gL, EL, Cm
Gate functions : alpha_n/m/h, beta_n/m/h
ODE : hh_rhs
Utilities : gates_ss, rest_ic, simulate_hh, find_limit_cycle_ic
"""
import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks

# ── Parameters ────────────────────────────────────────────────
gK  = 36.0;  EK  = -77.0
gNa = 120.0; ENa =  50.0
gL  =  0.3;  EL  = -54.4
Cm  =  1.0

# ── Gate rate functions ───────────────────────────────────────
def alpha_n(v): return 0.01*(v+55) / (1.0 - np.exp(-(v+55)/10.0))
def beta_n(v):  return 0.125 * np.exp(-(v+65)/80.0)
def alpha_m(v): return 0.1*(v+40)  / (1.0 - np.exp(-(v+40)/10.0))
def beta_m(v):  return 4.0  * np.exp(-(v+65)/18.0)
def alpha_h(v): return 0.07 * np.exp(-(v+65)/20.0)
def beta_h(v):  return 1.0  / (1.0 + np.exp(-(v+35)/10.0))

# ── ODE ───────────────────────────────────────────────────────
def hh_rhs(y, t, I_app):
    V, m, h, n = y
    dV = (I_app - gNa*m**3*h*(V-ENa) - gK*n**4*(V-EK) - gL*(V-EL)) / Cm
    return [dV,
            alpha_m(V)*(1-m) - beta_m(V)*m,
            alpha_h(V)*(1-h) - beta_h(V)*h,
            alpha_n(V)*(1-n) - beta_n(V)*n]

# ── Initial conditions ────────────────────────────────────────
def gates_ss(V):
    """Steady-state gate values alpha/(alpha+beta) at voltage V."""
    return (alpha_m(V)/(alpha_m(V)+beta_m(V)),
            alpha_h(V)/(alpha_h(V)+beta_h(V)),
            alpha_n(V)/(alpha_n(V)+beta_n(V)))

def rest_ic(V0=-65.):
    m0, h0, n0 = gates_ss(V0)
    return [V0, m0, h0, n0]

# ── Simulation ────────────────────────────────────────────────
def simulate_hh(I_app, T=100., dt=0.05, y0=None):
    """
    Simulate HH. Returns (t, V, m, h, n) as 1-D arrays.
    """
    t = np.arange(0, T+dt, dt)
    if y0 is None:
        y0 = rest_ic()
    sol = odeint(hh_rhs, y0, t, args=(I_app,))
    return t, sol[:,0], sol[:,1], sol[:,2], sol[:,3]

def find_limit_cycle_ic(I_app, gate_var="n", T_warmup=100., dt=0.05, n_skip=3):
    """
    Run HH from rest, skip the initial transient, and return a state
    already on the limit cycle.

    gate_var : "n" (xK) or "V" (xNa).
        For xK (slow gate n): trough is located on xK=n^4, because n
        lags V and anchoring on V's trough leaves xK mid-cycle.
        For xNa (fast gate m): trough is located on V directly, because
        m is so fast it tracks V nearly instantaneously.
    """
    t, V, m, h, n = simulate_hh(I_app, T_warmup, dt, y0=rest_ic())
    d = max(1, int(10/dt))
    peaks, _ = find_peaks(V, height=-20, distance=d)
    if len(peaks) < n_skip + 1:
        raise RuntimeError(
            f"only {len(peaks)} spikes in {T_warmup}ms warm-up; "
            f"need >= {n_skip+1} -- increase T_warmup")
    i_a, i_b = peaks[n_skip-1], peaks[n_skip]
    if gate_var == "n":
        xK = n**4
        i_trough = i_a + int(np.argmin(xK[i_a:i_b]))
    else:
        i_trough = i_a + int(np.argmin(V[i_a:i_b]))
    return [V[i_trough], m[i_trough], h[i_trough], n[i_trough]]
