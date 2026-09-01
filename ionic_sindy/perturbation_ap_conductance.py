"""
perturbation_apd_dvdt_conductance.py
=====================================
Tests the robustness of the HH and SINDy-hybrid models when the maximal
conductances gK and gNa are perturbed by ±20% from their nominal values.

Metrics evaluated on the LAST action potential out of 10 pulses:
    1. APD90  : action-potential duration at 90% repolarization (ms)
    2. dV/dt_max : maximum upstroke velocity (mV/ms)

Outputs (saved to outputs_apd_dvdt/):
    - traces_gK.png / traces_gNa.png     : all voltage traces (HH solid, SINDy dashed)
    - metrics_gK.png / metrics_gNa.png   : APD90 and dV/dt_max vs perturbation
    - error_gK.png   / error_gNa.png     : absolute and relative errors vs perturbation
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.integrate import odeint

import hh_model_v2 as hh


import scienceplots      
plt.style.use(["science", "ieee"])
plt.rc('font', size=14)

OUTDIR = "outputs_apd_dvdt"
os.makedirs(OUTDIR, exist_ok=True)

# ── Nominal conductances ──────────────────────────────────────
GK_NOM  = hh.gK    # 36.0
GNA_NOM = hh.gNa   # 120.0

# ── Perturbation levels ───────────────────────────────────────
PERTURBATIONS = np.array([-0.10, -0.08, -0.06, -0.04, -0.02, 0.0, +0.02, +0.04, +0.06, +0.08, +0.10])
LABELS  = [f"{int(round(p*100)):+d}%" for p in PERTURBATIONS]
NOM_IDX = int(np.argmin(np.abs(PERTURBATIONS)))

# ── Simulation settings ───────────────────────────────────────
# Run long enough to get at least 10 full pulses at I_app=10.
# Period ~14 ms at nominal => 10 pulses need ~140 ms + transient.
T       = 300.0   # ms -- more than enough
DT      = 0.01    # ms
I_TRAJ  = 12.0    # uA/cm2


# ── Color scheme (same as F-I curve scripts) ──────────────────
def _color(p):
    """Blue (negative) → black (nominal) → red (positive)."""
    if abs(p) < 1e-9:
        return (0.5, 0.5, 0.5)
    elif p < 0:
        t = abs(p) / 0.10
        return tuple((1-t)*np.array([0.75,0.87,1.0]) + t*np.array([0.05,0.27,0.55]))
    else:
        t = p / 0.10
        return tuple((1-t)*np.array([1.0,0.80,0.75]) + t*np.array([0.60,0.05,0.05]))


# =================================================================
# SINDy surrogate equations (fixed, identified at nominal g)
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
# Simulation functions
# =================================================================

def simulate_hh_g(I_app, gNa, gK, T=T, dt=DT, y0=None):
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
    sol = odeint(lambda y, t_: rhs(y, t_), y0, t)
    return t, sol[:, 0]


def simulate_sindy_g(I_app, gNa, gK, T=T, dt=DT, y0=None):
    if y0 is None:
        y0 = hh.rest_ic()
    V0, m0, h0, n0 = y0
    t = np.arange(0, T + dt, dt)
    V, xNa, xK = V0, m0**3 * h0, n0**4
    V_arr = np.zeros(len(t))
    V_arr[0] = V
    for i in range(1, len(t)):
        INa = gNa * xNa * (V - hh.ENa)
        IK  = gK  * xK  * (V - hh.EK)
        IL  = hh.gL     * (V - hh.EL)
        V   += dt * (I_app - INa - IK - IL) / hh.Cm
        xNa  = float(np.clip(xNa + dt * dxNa_dt(V, xNa), 0., 1.3))
        xK   = float(np.clip(xK  + dt * dxK_dt(V, xK),   0., 1.3))
        V_arr[i] = V
    return t, V_arr


# =================================================================
# Detect spike peaks and extract the last pulse window
# =================================================================

def detect_peaks(t, V, threshold=0.0, min_dist_ms=5.0):
    """Return indices of upward threshold crossings (one per AP)."""
    dt = t[1] - t[0]
    min_dist = max(1, int(min_dist_ms / dt))
    crossings = np.where((V[:-1] < threshold) & (V[1:] >= threshold))[0]
    # Remove crossings too close together
    if len(crossings) == 0:
        return np.array([], dtype=int)
    filtered = [crossings[0]]
    for c in crossings[1:]:
        if c - filtered[-1] >= min_dist:
            filtered.append(c)
    return np.array(filtered)


def extract_last_pulse_window(t, V, n_pulses=10, threshold=0.0,
                               discard_ms=100.0):
    """
    Find the LAST complete action potential in the steady-state portion
    of the trace (after discard_ms) and return a window covering that
    single pulse -- from its upward threshold crossing to just before
    the NEXT crossing (or to the end of the trace if no next crossing).

    Using the last spike in a time window (rather than counting exactly
    N crossings) is robust to phase drift between HH and SINDy: even
    when the two models are out of sync, each trace is evaluated on its
    OWN last complete spike.
    """
    mask = t >= discard_ms
    t_ss = t[mask]
    V_ss = V[mask]

    if len(t_ss) < 3:
        return None

    crossings = detect_peaks(t_ss, V_ss, threshold=threshold)
    if len(crossings) < 2:
        return None

    # Take the SECOND-TO-LAST crossing: this guarantees the window
    # ends at the NEXT (last) crossing, so the full repolarization
    # and AHP are always included.
    start_idx = crossings[-2]
    end_idx   = crossings[-1]          # end at next upward crossing

    t_win = t_ss[start_idx:end_idx]
    V_win = V_ss[start_idx:end_idx]

    if len(V_win) < 10 or np.max(V_win) < threshold:
        return None

    return t_win, V_win


# =================================================================
# Metrics: spike amplitude and V_AHP -- both threshold-free
# =================================================================

def compute_spike_amplitude(t_win, V_win):
    """
    Spike amplitude = V_peak - V_pre, where V_pre is the minimum
    voltage BEFORE the peak (the pre-spike resting level of this
    pulse). Threshold-free and robust to shifts in resting potential.
    """
    i_peak = int(np.argmax(V_win))
    V_peak = float(V_win[i_peak])
    V_pre  = float(np.min(V_win[:i_peak + 1]))
    return V_peak - V_pre, V_peak, V_pre


def compute_vahp(t_win, V_win):
    """
    Afterhyperpolarization: minimum voltage AFTER the peak.
    Sensitive to gK (drives repolarisation and AHP depth).
    """
    i_peak = int(np.argmax(V_win))
    if i_peak >= len(V_win) - 1:
        return np.nan
    return float(np.min(V_win[i_peak:]))


def compute_dvdt_max(t_win, V_win):
    """Maximum dV/dt (mV/ms) — maximum upstroke velocity."""
    dVdt = np.gradient(V_win, t_win)
    return float(np.max(dVdt))


# =================================================================
# Sweep
# =================================================================

def run_sweep(sweep_param, n_pulses=10):
    assert sweep_param in ("gK", "gNa")
    ic = hh.rest_ic()
    n_pert = len(PERTURBATIONS)

    amp_hh    = np.full(n_pert, np.nan)
    amp_sindy = np.full(n_pert, np.nan)
    vahp_hh   = np.full(n_pert, np.nan)
    vahp_sindy= np.full(n_pert, np.nan)
    dvdt_hh   = np.full(n_pert, np.nan)
    dvdt_sindy= np.full(n_pert, np.nan)

    traces_hh    = []
    traces_sindy = []
    t_ref = None

    for pi, p in enumerate(PERTURBATIONS):
        gK_use  = GK_NOM  * (1 + p) if sweep_param == "gK"  else GK_NOM
        gNa_use = GNA_NOM * (1 + p) if sweep_param == "gNa" else GNA_NOM

        t_hh, V_hh     = simulate_hh_g(I_TRAJ, gNa_use, gK_use, y0=ic)
        t_s,  V_sindy_ = simulate_sindy_g(I_TRAJ, gNa_use, gK_use, y0=ic)
        if t_ref is None:
            t_ref = t_hh

        traces_hh.append(V_hh)
        traces_sindy.append(V_sindy_)

        win_hh = extract_last_pulse_window(t_hh, V_hh, n_pulses=n_pulses)
        win_s  = extract_last_pulse_window(t_s, V_sindy_, n_pulses=n_pulses)

        if win_hh is not None:
            t_w, V_w = win_hh
            amp_hh[pi], _, _  = compute_spike_amplitude(t_w, V_w)
            vahp_hh[pi]       = compute_vahp(t_w, V_w)
            dvdt_hh[pi]       = compute_dvdt_max(t_w, V_w)
        if win_s is not None:
            t_w, V_w = win_s
            amp_sindy[pi], _, _ = compute_spike_amplitude(t_w, V_w)
            vahp_sindy[pi]      = compute_vahp(t_w, V_w)
            dvdt_sindy[pi]      = compute_dvdt_max(t_w, V_w)

        print(f"  {sweep_param} {LABELS[pi]:>4} | "
              f"HH:    amp={amp_hh[pi]:7.2f} mV, V_AHP={vahp_hh[pi]:7.2f} mV, "
              f"dVdt={dvdt_hh[pi]:7.2f} mV/ms | "
              f"SINDy: amp={amp_sindy[pi]:7.2f} mV, V_AHP={vahp_sindy[pi]:7.2f} mV, "
              f"dVdt={dvdt_sindy[pi]:7.2f} mV/ms")

    return (t_ref, traces_hh, traces_sindy,
            amp_hh, amp_sindy, vahp_hh, vahp_sindy,
            dvdt_hh, dvdt_sindy)


# =================================================================
# Plot 1: all voltage traces (HH solid, SINDy dashed)
# =================================================================

def plot_traces(t, traces_hh, traces_sindy, sweep_param, t_xlim=(50, 200)):
    fig, axes = plt.subplots(len(PERTURBATIONS), 1,
                              figsize=(8, 1.6 * len(PERTURBATIONS)),
                              sharex=True)
    mask = (t >= t_xlim[0]) & (t <= t_xlim[1])
    t_win = t[mask]

    for pi, ax in enumerate(axes):
        p = PERTURBATIONS[pi]
        c = _color(p)
        lw = 1.4
        ax.plot(t_win, traces_hh[pi][mask],    "-",  color=c, lw=lw, label="HH")
        ax.plot(t_win, traces_sindy[pi][mask], "--",  color=c, lw=lw, label="SINDy")
        ax.set_ylabel(f"{LABELS[pi]}", fontsize=8, rotation=0, labelpad=28)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2, lw=0.4)
        if pi == 0:
            ax.legend(fontsize=7, loc="upper right", ncol=2)

    param_label = r"$g_K$" if sweep_param == "gK" else r"$g_{Na}$"
    axes[0].set_title(f"Voltage traces — {param_label} perturbation", fontsize=9)
    axes[-1].set_xlabel("t (ms)", fontsize=9)

    plt.tight_layout()
    fname = os.path.join(OUTDIR, f"traces_{sweep_param}.pdf")
    plt.savefig(fname, dpi=300)
    plt.close(fig)
    print(f"  Saved: {fname}")


# =================================================================
# Combined figures: 2 rows (gK top, gNa bottom) x 3 cols
# (dV/dt_max | amplitude | V_AHP)
# One figure for quantities, one for relative errors (%)
# =================================================================

def _fill_panel(ax, pct, y_hh, y_sindy):
    """Plot HH (circles) and SINDy (squares) points in one panel."""
    for pi, p in enumerate(PERTURBATIONS):
        c  = _color(p)
        ms = 7 if pi == NOM_IDX else 5
        ms = 5
        if not np.isnan(y_hh[pi]):
            ax.plot(pct[pi], y_hh[pi],    "o", color=c, ms=ms, alpha=0.9)
        if not np.isnan(y_sindy[pi]):
            ax.plot(pct[pi], y_sindy[pi], "s", color=c, ms=ms, alpha=0.9)
    # ax.axvline(0, color="0.6", lw=0.8, ls=":")
    ax.set_xticks(np.arange(-10, 11, 5))
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.25)


def plot_combined_quantities(results_gK, results_gNa):
    """
    2 rows x 3 cols:
      row 0 = gK perturbation
      row 1 = gNa perturbation
      col 0 = dV/dt_max
      col 1 = spike amplitude
      col 2 = V_AHP
    """
    pct = PERTURBATIONS * 100
    (_, _, _, amp_hh_gK, amp_sindy_gK, vahp_hh_gK, vahp_sindy_gK,
     dvdt_hh_gK, dvdt_sindy_gK) = results_gK
    (_, _, _, amp_hh_gNa, amp_sindy_gNa, vahp_hh_gNa, vahp_sindy_gNa,
     dvdt_hh_gNa, dvdt_sindy_gNa) = results_gNa

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)

    panels = [
        # row 0: gK
        (axes[0,0], dvdt_hh_gK,  dvdt_sindy_gK,  r"$\max(dV/dt)$ (mV/ms)", r"$g_K$ perturbation"),
        (axes[0,1], amp_hh_gK,   amp_sindy_gK,   "Spike amplitude (mV)",    r"$g_K$ perturbation"),
        (axes[0,2], vahp_hh_gK,  vahp_sindy_gK,  r"$V_{AHP}$ (mV)",         r"$g_K$ perturbation"),
        # row 1: gNa
        (axes[1,0], dvdt_hh_gNa, dvdt_sindy_gNa, r"$\max(dV/dt)$ (mV/ms)", r"$g_{Na}$ perturbation"),
        (axes[1,1], amp_hh_gNa,  amp_sindy_gNa,  "Spike amplitude (mV)",    r"$g_{Na}$ perturbation"),
        (axes[1,2], vahp_hh_gNa, vahp_sindy_gNa, r"$V_{AHP}$ (mV)",         r"$g_{Na}$ perturbation"),
    ]

    for ax, y_hh, y_sindy, ylabel, row_label in panels:
        _fill_panel(ax, pct, y_hh, y_sindy)
        ax.set_ylabel(ylabel, fontsize=9)

    # Row labels on left column
    # for row, label in enumerate([r"$g_K$ perturbation", r"$g_{Na}$ perturbation"]):
    #     axes[row, 0].annotate(label, xy=(-0.35, 0.5), xycoords="axes fraction",
    #                            fontsize=9, ha="center", va="center", rotation=90)

    axes[0, 1].set_title(r"$g_K$ perturbation",  loc="center", fontsize=10, fontweight="bold")
    axes[1, 1].set_title(r"$g_{Na}$ perturbation", loc="center", fontsize=10, fontweight="bold")

    # x-axis labels on bottom row only
    for ax in axes[1, :]:
        ax.set_xlabel("Perturbation ($\%$)", fontsize=9)

    # # Column titles on top row only
    # for ax, title in zip(axes[0, :], [r"$\max(dV/dt)$", "Spike amplitude", r"$V_{AHP}$"]):
    #     ax.set_title(title, fontsize=10)

    # Legend once, top-right panel
    hh_h    = mlines.Line2D([], [], color="k", marker="o", ls="None", ms=6, label="HH")
    sindy_h = mlines.Line2D([], [], color="k", marker="s", ls="None", ms=6, label="SINDy")
    axes[0, 2].legend(handles=[hh_h, sindy_h], fontsize=8, framealpha=0.85,
                      loc="upper right")

    plt.tight_layout()
    fname = os.path.join(OUTDIR, "combined_quantities.pdf")
    plt.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fname}")


def plot_combined_errors(results_gK, results_gNa):
    """
    Same 2x3 layout but showing relative error (%) between SINDy and HH.
    """
    pct = PERTURBATIONS * 100
    (_, _, _, amp_hh_gK, amp_sindy_gK, vahp_hh_gK, vahp_sindy_gK,
     dvdt_hh_gK, dvdt_sindy_gK) = results_gK
    (_, _, _, amp_hh_gNa, amp_sindy_gNa, vahp_hh_gNa, vahp_sindy_gNa,
     dvdt_hh_gNa, dvdt_sindy_gNa) = results_gNa

    def rel_err(a, b):
        e = np.abs(a - b)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(np.abs(b) > 1e-9, 100.0 * e / np.abs(b), np.nan)
        return r

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)

    panels = [
        (axes[0,0], rel_err(dvdt_sindy_gK,  dvdt_hh_gK),
         r"$\max(dV/dt)$ error (%)", r"$g_K$ perturbation"),
        (axes[0,1], rel_err(amp_sindy_gK,   amp_hh_gK),
         "Amplitude error (%)",       r"$g_K$ perturbation"),
        (axes[0,2], rel_err(vahp_sindy_gK,  vahp_hh_gK),
         r"$V_{AHP}$ error (%)",      r"$g_K$ perturbation"),
        (axes[1,0], rel_err(dvdt_sindy_gNa, dvdt_hh_gNa),
         r"$\max(dV/dt)$ error (%)", r"$g_{Na}$ perturbation"),
        (axes[1,1], rel_err(amp_sindy_gNa,  amp_hh_gNa),
         "Amplitude error (%)",       r"$g_{Na}$ perturbation"),
        (axes[1,2], rel_err(vahp_sindy_gNa, vahp_hh_gNa),
         r"$V_{AHP}$ error (%)",      r"$g_{Na}$ perturbation"),
    ]

    for ax, err, ylabel, _ in panels:
        for pi, p in enumerate(PERTURBATIONS):
            if not np.isnan(err[pi]):
                c  = _color(p)
                ms = 7 if pi == NOM_IDX else 5
                ax.plot(pct[pi], err[pi], "D", color=c, ms=ms, alpha=0.9)
        # ax.axvline(0, color="0.6", lw=0.8, ls=":")
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xticks(np.arange(-10, 11, 5))
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.25)

    # for row, label in enumerate([r"$g_K$ perturbation", r"$g_{Na}$ perturbation"]):
        # axes[row, 0].annotate(label, xy=(-0.35, 0.5), xycoords="axes fraction",
                            #    fontsize=9, ha="center", va="center", rotation=90)

    axes[0, 1].set_title(r"$g_K$ perturbation",  loc="center", fontsize=10, fontweight="bold")
    axes[1, 1].set_title(r"$g_{Na}$ perturbation", loc="center", fontsize=10, fontweight="bold")

    for ax in axes[1, :]:
        ax.set_xlabel("Perturbation ($\%$)", fontsize=9)

    # for ax, title in zip(axes[0, :], [r"$\max(dV/dt)$", "Spike amplitude", r"$V_{AHP}$"]):
        # ax.set_title(title, fontsize=10)

    fig.suptitle("|SINDy − HH| relative error (%)", fontsize=11)
    plt.tight_layout()
    fname = os.path.join(OUTDIR, "combined_errors.pdf")
    plt.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fname}")


# =================================================================
# Main
# =================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Spike amplitude, V_AHP, dV/dt_max — using LAST of 10 pulses")
    print("=" * 70)

    results = {}
    for sweep_param in ("gK", "gNa"):
        print(f"\n--- Sweeping {sweep_param} ---")
        res = run_sweep(sweep_param, n_pulses=10)
        results[sweep_param] = res
        t_ref, traces_hh, traces_sindy = res[0], res[1], res[2]
        plot_traces(t_ref, traces_hh, traces_sindy, sweep_param)

    print("\n--- Generating combined figures ---")
    plot_combined_quantities(results["gK"], results["gNa"])
    plot_combined_errors(results["gK"], results["gNa"])

    print("\nDone.")
