"""
SINDy for xNa = m^3*h (sodium gating, Hodgkin-Huxley).

Observable:  xNa = INa / (gNa*(V-ENa)) = m^3*h
Derivative:  d(m^3*h)/dt = 3*m^2*(dm/dt)*h + m^3*(dh/dt)  -- analytic
Library:     Heaviside-gated: H(V_gate-V)*V^k*xNa, H(V_gate-V)*V^k*(1-xNa),
                               H(V-V_gate)*V^k*xNa, H(V-V_gate)*V^k*(1-xNa)
             V_gate = -40 mV  (m_inf=0.5, best from sweep)
             Uses IdentityLibrary + pre-computed Theta (make_theta).
Optimizer:   ps.STLSQ or ps.SR3 (from sindy_tools)
"""
import os, warnings
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots          # noqa: F401
plt.style.use(["science", "ieee", "no-latex"])
import pysindy as ps
from scipy.signal import find_peaks
from sklearn.metrics import r2_score, mean_squared_error
import re

from hh_model_v2 import (gK, EK, gNa, ENa, gL, EL, Cm,
                          alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h,
                          simulate_hh, rest_ic, find_limit_cycle_ic)
from sindy_tools import format_equation

os.makedirs("outputs_hh_xna_hs", exist_ok=True)

# ── Observable and analytic derivative ────────────────────────
def compute_x(m, h, n):
    return m**3 * h

def compute_dx(V, m, h, n):
    dm = alpha_m(V)*(1-m) - beta_m(V)*m
    dh = alpha_h(V)*(1-h) - beta_h(V)*h
    return 3.*m**2*dm*h + m**3*dh

# ── Heaviside library ─────────────────────────────────────────
def H(x):
    return np.heaviside(x, 0.5)

def make_theta(V, xNa, degree, V_gate):
    """
    Builds Theta with BOTH Heaviside halves:
      H(V_gate-V)*V^k*xNa, H(V_gate-V)*V^k*(1-xNa)  [below: recovery]
      H(V-V_gate)*V^k*xNa, H(V-V_gate)*V^k*(1-xNa)  [above: activation]
    Total: 4*(degree+1) columns.
    """
    Hm = H(V_gate - V)
    Hp = H(V - V_gate)
    cols, names = [], []
    for Hv, tag in [(Hm, f"H({V_gate:.0f}-V)"), (Hp, f"H(V-{V_gate:.0f})")]:
        for k in range(degree+1):
            Vk = V**k
            cols.append(Hv * Vk * xNa)
            cols.append(Hv * Vk * (1 - xNa))
            pfx = tag if k == 0 else f"{tag}*V^{k}"
            names.append(f"{pfx}*xNa")
            names.append(f"{pfx}*(1-xNa)")
    return np.column_stack(cols), names

# ── Fit, evaluate, R2  ────────────────────────────────────────
def fit_sindy(V, xNa, dxNa, degree, V_gate, threshold, opt_name="STLSQ"):
    Theta, feat_names = make_theta(V, xNa, degree, V_gate)
    Xdot = np.column_stack([np.zeros_like(dxNa), dxNa])

    if opt_name == "STLSQ":
        opt = ps.STLSQ(threshold=threshold, max_iter=2000, normalize_columns=True)
    else:
        opt = ps.SR3(reg_weight_lam=threshold, regularizer="l0",
                     max_iter=2000, normalize_columns=True)

    model = ps.SINDy(feature_library=ps.IdentityLibrary(), optimizer=opt)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        model.fit(Theta, t=1.0, x_dot=Xdot)

    c = model.coefficients()[1]
    if np.all(np.abs(c) < 1e-10):
        return None, feat_names
    return model, feat_names

def eval_dxNa(model, V_gate, degree, V_scalar, xNa_scalar):
    Th, _ = make_theta(np.array([V_scalar]), np.array([xNa_scalar]),
                       degree, V_gate)
    return float(Th[0] @ model.coefficients()[1])

def compute_r2(model, V_gate, degree, V, xNa, dxNa):
    Th, _ = make_theta(V, xNa, degree, V_gate)
    c = model.coefficients()[1]
    return r2_score(dxNa, Th @ c), int(np.sum(np.abs(c) > 1e-10))

# ── Compact equation print ─────────────────────────────────────
def print_compact_equation(model, feat_names, V_gate):
    """
    Prints the identified equation in factored alpha/beta form:
        dxNa/dt = H(V_gate-V)*[alpha_d(V)*(1-xNa) - beta_d(V)*xNa]
                + H(V-V_gate)*[alpha_p(V)*(1-xNa) - beta_p(V)*xNa]
    """
    c = model.coefficients()[1]
    tag_m = f"H({V_gate:.0f}-V)"
    tag_p = f"H(V-{V_gate:.0f})"
    buckets = {"Hm_xna": {}, "Hm_1mxna": {}, "Hp_xna": {}, "Hp_1mxna": {}}

    for coef, name in zip(c, feat_names):
        if abs(coef) < 1e-10:
            continue
        if name.startswith(tag_m):
            regime, rest = "Hm", name[len(tag_m):]
        elif name.startswith(tag_p):
            regime, rest = "Hp", name[len(tag_p):]
        else:
            print(f"  WARNING: cannot parse term '{name}'."); return
        if rest in ("*xNa", "xNa"):
            power, gating = 0, "xna"
        elif rest in ("*(1-xNa)", "(1-xNa)"):
            power, gating = 0, "1mxna"
        else:
            m_ = re.match(r"\*V\^(\d+)\*(.*)", rest)
            if m_ is None:
                print(f"  WARNING: cannot parse term '{name}'."); return
            power = int(m_.group(1))
            gating = "xna" if m_.group(2) == "xNa" else "1mxna"
        buckets[f"{regime}_{gating}"][power] = coef

    def poly_str(coeffs):
        if not coeffs: return "0"
        terms = []
        for k in sorted(coeffs):
            c_ = coeffs[k]
            if   k == 0: terms.append(f"{c_:+.4e}")
            elif k == 1: terms.append(f"{c_:+.4e}*V")
            else:        terms.append(f"{c_:+.4e}*V^{k}")
        return "  ".join(terms)

    bd = {k: -c_ for k, c_ in buckets["Hm_xna"].items()}
    bp = {k: -c_ for k, c_ in buckets["Hp_xna"].items()}
    print(f"\n{'='*60}")
    print("COMPACT FORM (alpha/beta):")
    print(f"  dxNa/dt = H({V_gate:.0f}-V)*[alpha_d(V)*(1-xNa) - beta_d(V)*xNa]")
    print(f"          + H(V-{V_gate:.0f})*[alpha_p(V)*(1-xNa) - beta_p(V)*xNa]")
    print(f"\n  alpha_d(V) = {poly_str(buckets['Hm_1mxna'])}")
    print(f"   beta_d(V) = {poly_str(bd)}")
    print(f"\n  alpha_p(V) = {poly_str(buckets['Hp_1mxna'])}")
    print(f"   beta_p(V) = {poly_str(bp)}")

# ── Simulate identified system ─────────────────────────────────
def simulate_id(model, degree, V_gate, I_app=10., T=100., dt=0.05, y0=None):
    """
    HH hybrid: V and n classical. xNa from SINDy.
    INa = gNa * xNa * (V-ENa)  [xNa from SINDy]
    IK  = gK  * n^4 * (V-EK)   [classical]
    """
    t = np.arange(0, T+dt, dt)
    if y0 is None:
        y0 = rest_ic()
    V0, m0, h0, n0 = y0
    xNa0 = m0**3 * h0

    V, n, xNa = V0, n0, xNa0
    V_arr   = np.zeros(len(t)); V_arr[0]   = V
    xNa_arr = np.zeros(len(t)); xNa_arr[0] = xNa

    for i in range(1, len(t)):
        INa = gNa * xNa  * (V - ENa)
        IK  = gK  * n**4 * (V - EK)
        IL  = gL          * (V - EL)
        dV  = (I_app - INa - IK - IL) / Cm
        dn  = alpha_n(V)*(1-n) - beta_n(V)*n
        dxNa_dt = eval_dxNa(model, V_gate, degree, V, xNa)

        V   += dt*dV; n += dt*dn
        xNa  = float(np.clip(xNa + dt*dxNa_dt, 0., 1.3))
        V_arr[i] = V; xNa_arr[i] = xNa

    return t, V_arr, xNa_arr

# ── Main ──────────────────────────────────────────────────────
if __name__ == "__main__":

    I_APP      = 10.0
    T          = 99.0
    DT         = 0.05
    DEGREES    = [2, 3, 4]
    V_GATE     = -40.0
    OPT_NAME   = "SR3"
    THRESHOLDS = [1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.3, 0.5, 1.0]

    print("HH SINDy xNa = m^3*h | Heaviside library\n")

    ic = find_limit_cycle_ic(I_APP, gate_var="V", T_warmup=100., dt=DT)
    print(f"limit-cycle IC: V0={ic[0]:.2f} m0={ic[1]:.4f} "
          f"h0={ic[2]:.4f} n0={ic[3]:.4f}\n")

    t, V, m, h, n = simulate_hh(I_APP, T, DT, y0=ic)
    xNa  = compute_x(m, h, n)
    dxNa = compute_dx(V, m, h, n)
    print(f"V   : [{V.min():.1f}, {V.max():.1f}] mV")
    print(f"xNa : [{xNa.min():.4f}, {xNa.max():.4f}]")
    print(f"max|dxNa/dt| = {np.abs(dxNa).max():.4f} ms-1\n")

    hdr = "deg | thresh   |   n | R2_fit   | RMSE_V   | picos"
    print(hdr); print("-"*len(hdr))

    best, best_rmse = None, np.inf
    results = []

    for degree in DEGREES:
        for th in THRESHOLDS:
            try:
                model, feat_names = fit_sindy(V, xNa, dxNa, degree, V_GATE,
                                               th, OPT_NAME)
                if model is None:
                    print(f"  {degree} | {th:7.0e} |   0 | -- nulo"); continue

                r2, n_nz = compute_r2(model, V_GATE, degree, V, xNa, dxNa)
                t_id, V_id, xNa_id = simulate_id(
                    model, degree, V_GATE, I_app=I_APP, T=T, dt=DT, y0=ic)

                rmse_v = float(np.sqrt(mean_squared_error(V, V_id))) \
                         if np.all(np.isfinite(V_id)) else np.nan
                d   = max(1, int(10/DT))
                n_o = len(find_peaks(V,    height=-20, distance=d)[0])
                n_i = len(find_peaks(V_id, height=-20, distance=d)[0])

                print(f"  {degree} | {th:7.0e} | {n_nz:3} |"
                      f" {r2:.5f} | {rmse_v:.4f} | {n_o}/{n_i}")

                if np.isfinite(rmse_v) and rmse_v < best_rmse:
                    best_rmse = rmse_v
                    best = dict(degree=degree, th=th,
                                model=model, feat_names=feat_names,
                                r2=r2, n_nz=n_nz,
                                t_id=t_id, V_id=V_id, xNa_id=xNa_id)
                results.append(dict(degree=degree, th=th,
                                    model=model, feat_names=feat_names,
                                    r2=r2, n_nz=n_nz, rmse_v=rmse_v,
                                    t_id=t_id, V_id=V_id, xNa_id=xNa_id))
            except Exception as e:
                print(f"  {degree} | {th:.0e} | ERRO: {e}")

    # --- equations ---
    print(f"\n{'='*65}\nEQUATIONS\n{'='*65}")
    for r in results:
        print(f"\n{'='*60}\nIDENTIFIED MODEL:")
        print(f"\ndeg={r['degree']} th={r['th']:.0e}"
              f"  R2={r['r2']:.4f}  RMSE_V={r['rmse_v']:.4f}  n={r['n_nz']}")
        print(f"  dxNa/dt = {format_equation(r['model'], r['feat_names'])}")
        print_compact_equation(r["model"], r["feat_names"], V_GATE)

    # --- plots ---
    colors_per = plt.cm.tab10(np.linspace(0, 1, 10))
    for degree in DEGREES:
        dr = [r for r in results if r["degree"] == degree]
        if not dr: continue

        # overview
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].plot(t, V,    "--k", lw=1.5, label="original", zorder=10)
        axes[1].plot(t, xNa, "--k", lw=1.5, label="original", zorder=10)
        colors = plt.cm.tab10(np.linspace(0, 1, len(dr)))
        for r, col in zip(dr, colors):
            lbl = f"th={r['th']:.0e} (n={r['n_nz']})"
            axes[0].plot(r["t_id"], r["V_id"],   color=col, lw=1.2, alpha=0.85, label=lbl)
            axes[1].plot(r["t_id"], r["xNa_id"], color=col, lw=1.2, alpha=0.85, label=lbl)
        axes[0].set(xlabel="t (ms)", ylabel="V (mV)",
                    title=f"V(t) deg={degree}")
        axes[1].set(xlabel="t (ms)", ylabel="xNa",
                    title=f"xNa(t) deg={degree}")
        axes[0].legend(fontsize=7); axes[1].legend(fontsize=7)
        plt.suptitle(f"HH SINDy xNa — Heaviside  V_gate={V_GATE:.0f}mV  deg={degree}",
                     fontsize=11)
        plt.tight_layout()
        fname = f"outputs_hh_xna_hs/xna_hs_deg{degree}.png"
        plt.savefig(fname, dpi=140); plt.close()
        print(f"  Saved: {fname}")

        # one PNG per model
        for r in dr:
            fig, axes = plt.subplots(1, 2, figsize=(7, 3))
            axes[0].plot(t, V,         "--k", lw=1.5, label="original")
            axes[0].plot(r["t_id"], r["V_id"],   color=colors_per[0], lw=1.2, label="SINDy")
            axes[1].plot(t, xNa,       "--k", lw=1.5, label="original")
            axes[1].plot(r["t_id"], r["xNa_id"], color=colors_per[2], lw=1.2, label="SINDy")
            axes[0].set(xlabel="t (ms)", ylabel="V (mV)")
            axes[1].set(xlabel="t (ms)", ylabel="xNa")
            axes[0].legend(fontsize=7, frameon=True, loc=1)
            axes[1].legend(fontsize=7, frameon=True, loc=1)
            plt.suptitle(
                f"xNa  deg={r['degree']}  th={r['th']:.0e}"
                f"  n={r['n_nz']}  R2={r['r2']:.4f}  RMSE={r['rmse_v']:.2f}mV",
                fontsize=8)
            plt.tight_layout()
            fname_r = (f"outputs_hh_xna_hs/"
                       f"hh_sindy_xna_deg{r['degree']}_th{r['th']:.0e}.png")
            plt.savefig(fname_r, dpi=140); plt.close()
            print(f"Saved: {fname_r}")

    if best:
        print(f"\nBEST: deg={best['degree']} th={best['th']:.0e}"
              f"  RMSE_V={best_rmse:.4f}  R2={best['r2']:.5f}  n={best['n_nz']}"
              f"  V_gate={V_GATE:.0f}mV")
        print("  dxNa/dt = " + format_equation(best["model"], best["feat_names"]))
        print_compact_equation(best["model"], best["feat_names"], V_GATE)


# --- Export metrics (R2, RMSE) vs threshold to .txt ---

    metrics_dir = "outputs_hh_xna_hs/metrics_txt"
    os.makedirs(metrics_dir, exist_ok=True)

    # Um arquivo por degree (útil se quiser comparar graus separadamente)
    for degree in DEGREES:
        dr = sorted([r for r in results if r["degree"] == degree],
                    key=lambda r: r["th"])
        if not dr:
            continue

        fname = f"{metrics_dir}/metrics_deg{degree}.txt"
        with open(fname, "w") as f:
            f.write("# threshold\tRMSE_V\tR2\tn_nonzero\n")
            for r in dr:
                f.write(f"{r['th']:.6e}\t{r['rmse_v']:.6e}\t"
                        f"{r['r2']:.6e}\t{r['n_nz']}\n")
        print(f"  Saved: {fname}")

    # Arquivo único com todos os degrees juntos (coluna extra "degree")
    fname_all = f"{metrics_dir}/metrics_all_degrees.txt"
    with open(fname_all, "w") as f:
        f.write("# degree\tthreshold\tRMSE_V\tR2\tn_nonzero\n")
        for r in sorted(results, key=lambda r: (r["degree"], r["th"])):
            f.write(f"{r['degree']}\t{r['th']:.6e}\t{r['rmse_v']:.6e}\t"
                    f"{r['r2']:.6e}\t{r['n_nz']}\n")
    print(f"  Saved: {fname_all}")
