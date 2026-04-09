import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# =============================================================================
# FITZHUGH-NAGUMO MODEL 
# =============================================================================

def fhn(t, z):
    x, y = z[0], z[1]
    dxdt = x - y - (1.0/3.0)*x*x*x + 0.5  # fast variable (membrane potential)
    dydt = 0.040*x - 0.028*y + 0.032      # slow variable (recovery)
    return (dxdt, dydt)

# Vectorized version (used for derivative computation) 
def fhn_vectorized(t, z):
    dydt = np.empty_like(z)
    x, y = z[0], z[1]
    dydt[0, ...] = x - y - (1.0/3.0)*x*x*x + 0.5  # dx/dt for all timesteps
    dydt[1, ...] = 0.040*x - 0.028*y + 0.032      # dy/dt for all timesteps
    return dydt

# =============================================================================
# SOLVER — Adaptive integration via odeint (LSODA)
# =============================================================================

def solve_fhn(y0, t_final=300.0, dt=0.01):

    # Time axis
    dt      = 0.01
    t_final = 300
    nsteps  = int(t_final / dt) + 1
    t       = np.linspace(0, t_final, nsteps)

    # Solver settings
    dt_max = t_final / 1000             # max internal step size
    y0     = (0.0, 0.0)                 # initial conditions: (x0, y0)
    rtol   = 10e-12                     # relative tolerance
    atol   = 10e-12 * np.ones_like(y0)  # absolute tolerance (per variable)

    # Integrate 
    y = odeint(fhn, y0, t, tfirst=True, rtol=rtol, atol=atol, hmax=dt_max)

    # Compute derivatives at each timestep
    dy = fhn_vectorized(t, y.T)  # expects (2, nsteps), returns (2, nsteps)
    dy = dy.T                    # back to (nsteps, 2)

    # Pack output as (5, N) array 
    # rows: t | x | y | dx/dt | dy/dt
    data = np.vstack((t, y[:, 0], y[:, 1], dy[:, 0], dy[:, 1]))
    return data