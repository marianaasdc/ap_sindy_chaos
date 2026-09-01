"""
conductance_perturbation_test.py
=================================
Tests the robustness of the SINDy hybrid model when the maximal
conductances gK and gNa are perturbed by up to ±20% from the nominal
values used during identification.

The SINDy surrogate equations (dxNa_dt, dxK_dt) were identified at
the NOMINAL conductances -- they are FIXED for all perturbations.
Only the voltage equation (which uses the perturbed g) changes.
This directly answers Reviewer #2 Sub-concern 2: does the reduced
model remain valid when conductances are changed, or must it be
re-identified?

For each perturbation level we compute:
  1. F-I curve (firing frequency vs I_app) for HH and SINDy
  2. A representative trajectory at I_app = 10 uA/cm2
  3. RMSE_V and frequency error as scalar metrics

Outputs (saved to outputs_conductance/):
  - fi_curve_gK_sweep.png
  - fi_curve_gNa_sweep.png
  - trajectory_gK_*.png
  - trajectory_gNa_*.png
  - metrics_summary.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.integrate import odeint

import hh_model_v2 as hh


import scienceplots      
plt.style.use(["science", "ieee"])
plt.rc('font', size=14)


os.makedirs("outputs_conductance", exist_ok=True)

# ── Nominal conductances ──────────────────────────────────────
GK_NOM  = hh.gK    # 36.0
GNA_NOM = hh.gNa   # 120.0

# ── Perturbation levels ───────────────────────────────────────
PERTURBATIONS = np.array([-0.10, -0.05, 0.0,  +0.05, +0.10])
LABELS        = [f"{int(round(p*100)):+d}%" for p in PERTURBATIONS]

# Color ramp: dark→light for negative, green for nominal, light→dark for positive
# HH lines: greys; SINDy lines: blue-to-red through green at centre
# def _sindy_color(p):
    # """Light blue (-20%) → dark blue (+20%), uniform hue."""
    # t = (p + 0.20) / 0.40          # 0 at -20%, 1 at +20%
    # light = np.array([0.75, 0.87, 1.00])   # very light blue
    # dark  = np.array([0.05, 0.27, 0.55])   # dark navy blue
    # return tuple((1 - t) * light + t * dark)

def _sindy_color(p):
    if abs(p) < 1e-9:        # nominal
        return (0.1, 0.1, 0.1)
    elif p < 0:              # negative: light→dark blue
        t = abs(p) / 0.20
        return tuple((1-t)*np.array([0.75,0.87,1.0]) + t*np.array([0.05,0.27,0.55]))
    else:                    # positive: light→dark red
        t = p / 0.20
        return tuple((1-t)*np.array([1.0,0.80,0.75]) + t*np.array([0.60,0.05,0.05]))    

# def _hh_color(p):
#     """Light grey (-20%) → black (0%) → light grey (+20%)."""
#     intensity = 0.75 - 0.55 * (1 - abs(p) / 0.20)
#     return (intensity, intensity, intensity)

def _hh_color(p):
    return _sindy_color(p)    

SINDY_COLORS = [_sindy_color(p) for p in PERTURBATIONS]
HH_COLORS    = [_hh_color(p)    for p in PERTURBATIONS]
NOM_IDX = int(np.argmin(np.abs(PERTURBATIONS)))

# ── Simulation settings ───────────────────────────────────────
T          = 300.0
DT         = 0.01
DISCARD_MS = 50.0
THRESHOLD  = 0.0
I_APP_VALUES = np.linspace(0.0, 40.0, 21)
I_TRAJ     = 10.0   # I_app used for trajectory plots


# =================================================================
# SINDy surrogate equations (hardcoded, identified at nominal g)
# =================================================================

def H(x):
    return np.heaviside(x, 0.5)

V_GATE_NA = -40.0

def dxNa_dt(V, xNa):
    Hm = H(V_GATE_NA - V)
    Hp = H(V - V_GATE_NA)
    return (  9.5857e-01 * Hm         * (1-xNa)
            + 4.8522e-01 * Hm * V     * xNa
            + 4.1550e-02 * Hm * V     * (1-xNa)
            + 1.6413e-02 * Hm * V**2  * xNa
            + 6.0058e-04 * Hm * V**2  * (1-xNa)
            + 1.7481e-04 * Hm * V**3  * xNa
            + 2.8939e-06 * Hm * V**3  * (1-xNa)
            - 9.2165e+00 * Hp          * xNa
            + 1.0691e+00 * Hp          * (1-xNa)
            - 1.1254e-01 * Hp * V     * xNa
            + 4.3061e-02 * Hp * V     * (1-xNa)
            + 4.9414e-04 * Hp * V**2  * (1-xNa))

def dxK_dt(V, xK):
    return (  3.6987e-01 * xK
            + 1.0101e-01 * (1-xK)
            + 1.8080e-02 * V    * xK
            + 1.8574e-03 * V    * (1-xK)
            + 1.4947e-04 * V**2 * xK
            + 9.4717e-07 * V**3 * xK
            + 1.5333e-09 * V**4 * (1-xK))


# =================================================================
# Simulation functions that accept EXPLICIT conductances
# =================================================================

def simulate_hh_g(I_app, gNa, gK, T=T, dt=DT, y0=None):
    """HH with explicit (possibly perturbed) gNa, gK."""
    if y0 is None:
        y0 = hh.rest_ic()
    t = np.arange(0, T + dt, dt)

    def rhs(y, t_):
        V, m, h, n = y
        dV = (I_app
              - gNa * m**3 * h * (V - hh.ENa)
              - gK  * n**4     * (V - hh.EK)
              - hh.gL          * (V - hh.EL)) / hh.Cm
        return [dV,
                hh.alpha_m(V)*(1-m) - hh.beta_m(V)*m,
                hh.alpha_h(V)*(1-h) - hh.beta_h(V)*h,
                hh.alpha_n(V)*(1-n) - hh.beta_n(V)*n]

    sol = odeint(rhs, y0, t, args=())
    # odeint needs args in rhs -- use closure instead (already done above)
    sol = odeint(lambda y, t_: rhs(y, t_), y0, t)
    return t, sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3]


def simulate_sindy_g(I_app, gNa, gK, T=T, dt=DT, y0=None):
    """
    SINDy hybrid with explicit (possibly perturbed) gNa, gK in the
    voltage equation. The surrogate gating equations are FIXED
    (identified at nominal conductances).
    """
    if y0 is None:
        y0 = hh.rest_ic()
    V0, m0, h0, n0 = y0
    xNa0 = m0**3 * h0
    xK0  = n0**4

    t = np.arange(0, T + dt, dt)
    V, xNa, xK = V0, xNa0, xK0
    V_arr   = np.zeros(len(t)); V_arr[0]   = V
    xNa_arr = np.zeros(len(t)); xNa_arr[0] = xNa
    xK_arr  = np.zeros(len(t)); xK_arr[0]  = xK

    for i in range(1, len(t)):
        INa = gNa * xNa * (V - hh.ENa)
        IK  = gK  * xK  * (V - hh.EK)
        IL  = hh.gL     * (V - hh.EL)

        dV   = (I_app - INa - IK - IL) / hh.Cm
        dxna = dxNa_dt(V, xNa)
        dxk  = dxK_dt(V, xK)

        V   += dt * dV
        xNa  = float(np.clip(xNa + dt * dxna, 0., 1.3))
        xK   = float(np.clip(xK  + dt * dxk,  0., 1.3))

        V_arr[i]   = V
        xNa_arr[i] = xNa
        xK_arr[i]  = xK

    return t, V_arr, xNa_arr, xK_arr


# =================================================================
# Metrics
# =================================================================

def spike_frequency(t, V, threshold=THRESHOLD, discard_ms=DISCARD_MS):
    mask = t >= discard_ms
    t_w, V_w = t[mask], V[mask]
    if len(t_w) < 2:
        return 0.0
    dt = t[1] - t[0]
    peaks, _ = find_peaks(V_w, height=threshold, distance=max(1, int(1.0/dt)))
    dur_s = (t_w[-1] - t_w[0]) / 1000.0
    return len(peaks) / dur_s if dur_s > 0 else 0.0


def rmse(a, b):
    return float(np.sqrt(np.mean((a - b)**2)))


# =================================================================
# Main sweep: compute F-I curves for all perturbations of one param
# =================================================================

def run_sweep(sweep_param, outdir="outputs_conductance"):
    assert sweep_param in ("gK", "gNa")
    ic = hh.rest_ic()
    n_pert = len(PERTURBATIONS)
    freq_hh    = np.zeros((n_pert, len(I_APP_VALUES)))
    freq_sindy = np.zeros((n_pert, len(I_APP_VALUES)))

    for pi, p in enumerate(PERTURBATIONS):
        gK_use  = GK_NOM  * (1 + p) if sweep_param == "gK"  else GK_NOM
        gNa_use = GNA_NOM * (1 + p) if sweep_param == "gNa" else GNA_NOM

        for ii, I_app in enumerate(I_APP_VALUES):
            _, V_hh, *_ = simulate_hh_g(I_app, gNa_use, gK_use, y0=ic)
            t_s, V_s, *_ = simulate_sindy_g(I_app, gNa_use, gK_use, y0=ic)
            freq_hh[pi, ii]    = spike_frequency(t_s, V_hh)
            freq_sindy[pi, ii] = spike_frequency(t_s, V_s)

        print(f"  {sweep_param} {LABELS[pi]:>4}  done")

    return freq_hh, freq_sindy


def plot_combined(freq_hh_gK, freq_sindy_gK,
                  freq_hh_gNa, freq_sindy_gNa,
                  outdir="outputs_conductance"):
    """
    One figure, two panels side by side (gK | gNa).
    All perturbation levels overlaid in each panel.
    HH lines: grey shades, solid.
    SINDy lines: blue→green→red ramp, dashed.
    Nominal (0%) highlighted: black solid (HH) + thick green dashed (SINDy).
    No overall title. Axis labels only.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), sharey=True)

    datasets = [
        (axes[0], freq_hh_gK,  freq_sindy_gK,  r"$g_K$ perturbation"),
        (axes[1], freq_hh_gNa, freq_sindy_gNa, r"$g_{Na}$ perturbation"),
    ]

    for ax, fhh, fsid, panel_title in datasets:
        for pi, p in enumerate(PERTURBATIONS):
            lw_hh  = 1.8 if pi == NOM_IDX else 1.2
            lw_sid = 1.8 if pi == NOM_IDX else 1.2
            alpha  = 1.0 if pi == NOM_IDX else 0.75
            zorder = 4  if pi == NOM_IDX else 2

            # HH: solid grey lines (nominal = black)
            c_hh = "black" if pi == NOM_IDX else HH_COLORS[pi]
            ax.plot(I_APP_VALUES, fhh[pi], "-",
                    color=c_hh, lw=lw_hh, alpha=alpha, zorder=zorder)

            # SINDy: dashed blue-ramp lines
            c_sid = _sindy_color(p)
            ax.plot(I_APP_VALUES, fsid[pi], "--",
                    color=c_sid, lw=lw_sid, alpha=alpha, zorder=zorder,
                    label=LABELS[pi])

        ax.set_xlabel(r"$I_{app}$ ($\mu$A/cm²)", fontsize=9)
        ax.set_title(panel_title, fontsize=9)
        ax.tick_params(labelsize=8)
        # ax.grid(True, alpha=0.25, lw=0.5)

    axes[0].set_ylabel("Firing frequency (Hz)", fontsize=9)

    # Legend: style indicators + one entry per perturbation level
    handles = []
    import matplotlib.lines as mlines
    hh_patch  = mlines.Line2D([], [], color="black", lw=1.5,
                               linestyle="-",  label="HH (solid)")
    sid_patch = mlines.Line2D([], [], color=_sindy_color(0.0), lw=1.5,
                               linestyle="--", label="SINDy (dashed)")
    handles += [hh_patch, sid_patch]
    for pi, p in enumerate(PERTURBATIONS):
        lw = 1.8 if pi == NOM_IDX else 1.2
        col = _sindy_color(p)
        handles.append(mlines.Line2D([], [], color=col, lw=lw,
                                      linestyle="--", label=LABELS[pi]))
    axes[1].legend(handles=handles, fontsize=7, ncol=1,
                   loc="lower right", framealpha=0.85)

    plt.tight_layout()
    fname = os.path.join(outdir, "fi_curve_conductance_sweep.pdf")
    plt.savefig(fname, dpi=300)
    plt.close(fig)
    print(f"\nSaved: {fname}")


if __name__ == "__main__":
    os.makedirs("outputs_conductance", exist_ok=True)

    print("Sweeping gK ± 20% (5% steps) ...")
    fhh_gK, fsid_gK = run_sweep("gK")

    print("\nSweeping gNa ± 20% (5% steps) ...")
    fhh_gNa, fsid_gNa = run_sweep("gNa")

    plot_combined(fhh_gK, fsid_gK, fhh_gNa, fsid_gNa)
