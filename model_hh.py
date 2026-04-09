import numpy as np

# =============================================================================
# Gatings - alpha & Beta rate functions
# =============================================================================

# --- Sodium activation (m) ---
def alpha_m(V):
    return 0.1 * (V + 40) / (1 - np.exp(-(V + 40) / 10))

def beta_m(V):
    return 4.0 * np.exp(-0.0556 * (V + 65))

# --- Sodium inactivation (h) ---
def alpha_h(V):
    return 0.07 * np.exp(-0.05 * (V + 65))

def beta_h(V):
    return 1 / (1 + np.exp(-(V + 35) / 10))

# --- Potassium activation (n) ---
def alpha_n(V):
    return 0.01 * (V + 55) / (1 - np.exp(-(V + 55) / 10))

def beta_n(V):
    return 0.125 * np.exp(-0.0125 * (V + 65))

# =============================================================================
# HODGKIN-HUXLEY MODEL
# =============================================================================

def hh_model(y, I):
    V, m, h, n = y

    C_m  = 1.0    # membrane capacitance  (µF/cm²)
    g_Na = 120.0  # max Na conductance    (mS/cm²)
    g_K  = 36.0   # max K  conductance    (mS/cm²)
    g_L  = 0.3    # leak conductance      (mS/cm²)
    E_Na = 50.0   # Na reversal potential (mV)
    E_K  = -77.0  # K  reversal potential (mV)
    E_L  = -54.4  # leak reversal         (mV)

    # ODEs
    dVdt = (I - g_Na * m**3 * h * (V - E_Na) - g_K * n**4 * (V - E_K) - g_L * (V - E_L)) / C_m
    dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h
    dndt = alpha_n(V) * (1 - n) - beta_n(V) * n

    return np.array([dVdt, dmdt, dhdt, dndt])

# =============================================================================
# SOLVER — Forward Euler integration
# =============================================================================

def solve_hh(y0, I=10.0, t_final=100.0, dt=0.05):

    n_steps = int(t_final / dt)
    t = np.arange(0, t_final, dt)

    V  = np.zeros(n_steps)
    m  = np.zeros(n_steps)
    h  = np.zeros(n_steps)
    n  = np.zeros(n_steps)
    dV = np.zeros(n_steps)
    dm = np.zeros(n_steps)
    dh = np.zeros(n_steps)
    dn = np.zeros(n_steps)

    # Initial conditions
    V[0], m[0], h[0], n[0] = y0
    dV[0], dm[0], dh[0], dn[0] = hh_model(y0, I)

    # Time integration
    for i in range(1, n_steps):
        dydt = hh_model([V[i-1], m[i-1], h[i-1], n[i-1]], I)
        V[i] = V[i-1] + dt * dydt[0]
        m[i] = m[i-1] + dt * dydt[1]
        h[i] = h[i-1] + dt * dydt[2]
        n[i] = n[i-1] + dt * dydt[3]
        dV[i], dm[i], dh[i], dn[i] = dydt

    # Pack output as (9, N) array
    # rows: t | V | m | h | n | dV | dm | dh | dn
    data = np.vstack((t, V, m, h, n, dV, dm, dh, dn))
    return data