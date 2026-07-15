"""
SINDy for xK = n^4 (potassium gating, Hodgkin-Huxley).

Variante de hh_sindy_xk_v2.py: em vez de varrer o degree do library,
o degree fica fixo (=3) e varre-se o OPTIMIZER (STLSQ, SR3, LASSO, SSR),
cada um percorrendo a mesma lista de thresholds.

Observable:  xK = IK / (gK*(V-EK)) = n^4
Derivative:  d(n^4)/dt = 4*n^3 * (alpha_n*(1-n) - beta_n*n)  -- analytic
Library:     GeneralizedLibrary: V^k*xK e V^k*(1-xK), k=0..degree (degree fixo=3)
             (from sindy_tools.make_library)
Optimizer:   ps.STLSQ / ps.SR3 / LASSO / ps.SSR (from sindy_tools.fit_sindy)
"""
import os, warnings
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots          # noqa: F401
plt.style.use(["science", "ieee", "no-latex"])
from scipy.signal import find_peaks
from sklearn.metrics import mean_squared_error
import re

from hh_model_v2 import (gK, EK, gNa, ENa, gL, EL, Cm,
                          alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h,
                          simulate_hh, rest_ic, find_limit_cycle_ic)
from sindy_tools import make_library, fit_sindy, eval_dx, compute_r2, format_equation

os.makedirs("outputs_hh_xk_opt", exist_ok=True)

# ── Observable and analytic derivative ────────────────────────
def compute_x(m, h, n):
    return n**4

def compute_dx(V, m, h, n):
    dn = alpha_n(V)*(1-n) - beta_n(V)*n
    return 4.0 * n**3 * dn

# ── Compact equation print ─────────────────────────────────────
def print_compact_equation_xk(model, feat_names):
    """
    Prints the identified equation in factored alpha/beta form:
        dxK/dt = alpha(V)*(1-xK) - beta(V)*xK
    Feature names from GeneralizedLibrary: '1 x1', 'x0 x1', etc.
    """
    c = model.coefficients()[1]
    buckets = {"xk": {}, "1mxk": {}}

    for coef, name in zip(c, feat_names):
        if abs(coef) < 1e-10:
            continue
        name = name.strip()
        if name.endswith("(1-x1)"):
            gating = "1mxk"; vpart = name[:-len("(1-x1)")].strip()
        elif name.endswith("x1"):
            gating = "xk";   vpart = name[:-len("x1")].strip()
        else:
            print(f"  WARNING: cannot parse term '{name}'."); return
        if vpart in ("", "1"):           power = 0
        elif vpart == "x0":              power = 1
        else:
            m_ = re.match(r"x0\^(\d+)", vpart)
            if m_: power = int(m_.group(1))
            else:  print(f"  WARNING: cannot parse V power in '{name}'."); return
        buckets[gating][power] = coef

    def poly_str(coeffs):
        if not coeffs: return "0"
        terms = []
        for k in sorted(coeffs):
            c_ = coeffs[k]
            if   k == 0: terms.append(f"{c_:+.4e}")
            elif k == 1: terms.append(f"{c_:+.4e}*V")
            else:        terms.append(f"{c_:+.4e}*V^{k}")
        return "  ".join(terms)

    alpha_coeffs = buckets["1mxk"]
    beta_coeffs  = {k: -c_ for k, c_ in buckets["xk"].items()}
    print(f"\n{'='*60}")
    print("COMPACT FORM (alpha/beta):")
    print(f"  dxK/dt = alpha(V)*(1-xK) - beta(V)*xK")
    print(f"\n  alpha(V) = {poly_str(alpha_coeffs)}")
    print(f"   beta(V) = {poly_str(beta_coeffs)}")

# ── Simulate identified system ─────────────────────────────────
def simulate_id(model, lib, I_app=10., T=100., dt=0.05, y0=None):
    """
    HH hybrid: V, m, h classical. xK from SINDy.
    IK = gK * xK * (V-EK)   [xK from SINDy]
    INa = gNa * m^3 * h * (V-ENa)  [classical]
    """
    t = np.arange(0, T+dt, dt)
    if y0 is None:
        y0 = rest_ic()
    V0, m0, h0, n0 = y0
    xK0 = n0**4

    V, m, h, xK = V0, m0, h0, xK0
    V_arr  = np.zeros(len(t)); V_arr[0]  = V
    xK_arr = np.zeros(len(t)); xK_arr[0] = xK

    for i in range(1, len(t)):
        INa = gNa * m**3 * h * (V - ENa)
        IK  = gK  * xK       * (V - EK)
        IL  = gL              * (V - EL)
        dV  = (I_app - INa - IK - IL) / Cm
        dm  = alpha_m(V)*(1-m) - beta_m(V)*m
        dh  = alpha_h(V)*(1-h) - beta_h(V)*h
        dxK_dt = eval_dx(model, lib, V, xK)

        V  += dt*dV; m += dt*dm; h += dt*dh
        xK  = float(np.clip(xK + dt*dxK_dt, 0., 1.3))
        V_arr[i] = V; xK_arr[i] = xK

    return t, V_arr, xK_arr

# ── Main ──────────────────────────────────────────────────────
if __name__ == "__main__":

    I_APP      = 10.0
    T          = 90.0
    DT         = 0.05
    DEGREE     = 3                                   # grau fixo (era varrido antes)
    OPT_NAMES  = ["STLSQ", "SR3", "LASSO", "SSR"]     # varredura de optimizers
    THRESHOLDS = [1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.3, 0.5, 0.75, 1.0, 2.0]

    print(f"HH SINDy xK = n^4 | GeneralizedLibrary | degree={DEGREE} (fixo)\n")

    ic = find_limit_cycle_ic(I_APP, gate_var="n", T_warmup=105., dt=DT)
    print(f"limit-cycle IC: V0={ic[0]:.2f} m0={ic[1]:.4f} "
          f"h0={ic[2]:.4f} n0={ic[3]:.4f}\n")

    t, V, m, h, n = simulate_hh(I_APP, T, DT, y0=ic)
    xK  = compute_x(m, h, n)
    dxK = compute_dx(V, m, h, n)
    print(f"V  : [{V.min():.1f}, {V.max():.1f}] mV")
    print(f"xK : [{xK.min():.4f}, {xK.max():.4f}]")
    print(f"max|dxK/dt| = {np.abs(dxK).max():.4f} ms-1\n")

    hdr = "opt      | thresh   |   n | R2_fit   | RMSE_V   | picos"
    print(hdr); print("-"*len(hdr))

    lib = make_library(DEGREE)   # library única, degree fixo

    best, best_rmse = None, np.inf
    results = []

    for opt_name in OPT_NAMES:
        for th in THRESHOLDS:
            try:
                model, lib_fit = fit_sindy(V, xK, dxK, lib, th, opt_name)
                if model is None:
                    print(f"  {opt_name:8s} | {th:7.0e} |   0 | -- nulo"); continue

                r2, n_nz = compute_r2(model, lib_fit, V, xK, dxK)
                t_id, V_id, xK_id = simulate_id(model, lib_fit,
                                                  I_app=I_APP, T=T, dt=DT, y0=ic)

                rmse_v = float(np.sqrt(mean_squared_error(V, V_id))) \
                         if np.all(np.isfinite(V_id)) else np.nan
                d   = max(1, int(10/DT))
                n_o = len(find_peaks(V,    height=-20, distance=d)[0])
                n_i = len(find_peaks(V_id, height=-20, distance=d)[0])

                print(f"  {opt_name:8s} | {th:7.0e} | {n_nz:3} | {r2:.5f} | "
                      f"{rmse_v:.4f} | {n_o}/{n_i}")

                if np.isfinite(rmse_v) and rmse_v < best_rmse:
                    best_rmse = rmse_v
                    best = dict(opt_name=opt_name, th=th, model=model, lib=lib_fit,
                                r2=r2, n_nz=n_nz,
                                t_id=t_id, V_id=V_id, xK_id=xK_id)
                results.append(dict(opt_name=opt_name, th=th, model=model, lib=lib_fit,
                                    r2=r2, n_nz=n_nz, rmse_v=rmse_v,
                                    t_id=t_id, V_id=V_id, xK_id=xK_id))
            except Exception as e:
                print(f"  {opt_name:8s} | {th:.0e} | ERRO: {e}")

    # --- equations ---
    print(f"\n{'='*60}\nEQUATIONS\n{'='*60}")
    for r in results:
        print(f"\n{'='*60}")
        print("IDENTIFIED MODEL:")
        print(f"\nopt={r['opt_name']} th={r['th']:.0e}"
              f"  R2={r['r2']:.4f}  RMSE_V={r['rmse_v']:.4f}  n={r['n_nz']}")
        feat_names = r["model"].get_feature_names()
        print(f"  dxK/dt = {format_equation(r['model'], feat_names)}")
        print_compact_equation_xk(r["model"], feat_names)

    # --- plots ---
    colors_per = plt.cm.tab10(np.linspace(0, 1, 10))
    for opt_name in OPT_NAMES:
        dr = [r for r in results if r["opt_name"] == opt_name]
        if not dr: continue

        # overview
        fig, axes = plt.subplots(1, 2, figsize=(7, 3))
        axes[0].plot(t, V,  "--k", lw=1.5, label="original", zorder=10)
        axes[1].plot(t, xK, "--k", lw=1.5, label="original", zorder=10)
        colors = plt.cm.tab10(np.linspace(0, 1, len(dr)))
        for r, col in zip(dr, colors):
            lbl = f"th={r['th']:.0e} (n={r['n_nz']})"
            axes[0].plot(r["t_id"], r["V_id"],  color=col, lw=1.2, alpha=0.85, label=lbl)
            axes[1].plot(r["t_id"], r["xK_id"], color=col, lw=1.2, alpha=0.85, label=lbl)
        axes[0].set(xlabel="t (ms)", ylabel="V (mV)")
        axes[1].set(xlabel="t (ms)", ylabel="xK")
        axes[1].legend(fontsize=7, loc="upper left",
                       bbox_to_anchor=(1.01, 1), borderaxespad=0.)
        plt.suptitle(f"HH SINDy xK — optimizer={opt_name} (deg={DEGREE})")
        plt.tight_layout(); plt.subplots_adjust(right=0.78)
        fname = f"outputs_hh_xk_opt/hh_sindy_xk_{opt_name}.png"
        plt.savefig(fname, dpi=140); plt.close()
        print(f"\nSaved: {fname}")

        # one PNG per model
        for r in dr:
            fig, axes = plt.subplots(1, 2, figsize=(7, 3))
            axes[0].plot(t, V,         "--k", lw=1.5, label="original")
            axes[0].plot(r["t_id"], r["V_id"],  color=colors_per[0], lw=1.2, label="SINDy")
            axes[1].plot(t, xK,        "--k", lw=1.5, label="original")
            axes[1].plot(r["t_id"], r["xK_id"], color=colors_per[3], lw=1.2, label="SINDy")
            axes[0].set(xlabel="t (ms)", ylabel="V (mV)")
            axes[1].set(xlabel="t (ms)", ylabel="xK")
            axes[0].legend(fontsize=7, frameon=True, loc=1)
            axes[1].legend(fontsize=7, frameon=True, loc=1)
            plt.suptitle(
                f"xK  opt={r['opt_name']}  th={r['th']:.0e}"
                f"  n={r['n_nz']}  R2={r['r2']:.4f}  RMSE={r['rmse_v']:.2f}mV",
                fontsize=8)
            plt.tight_layout()
            fname_r = (f"outputs_hh_xk_opt/"
                       f"hh_sindy_xk_{r['opt_name']}_th{r['th']:.0e}.png")
            plt.savefig(fname_r, dpi=140); plt.close()
            print(f"Saved: {fname_r}")

    if best:
        print(f"\nBEST: opt={best['opt_name']} th={best['th']:.0e}"
              f"  RMSE_V={best_rmse:.4f}  R2={best['r2']:.5f}  n={best['n_nz']}")
        feat_names = best["model"].get_feature_names()
        print("  dxK/dt = " + format_equation(best["model"], feat_names))
        print_compact_equation_xk(best["model"], feat_names)

# --- Exporta métricas (R2, RMSE) vs threshold para .txt ---

    metrics_dir = "outputs_hh_xk_opt/metrics_txt"
    os.makedirs(metrics_dir, exist_ok=True)

    # Um arquivo por optimizer (útil se quiser comparar optimizers separadamente)
    for opt_name in OPT_NAMES:
        dr = sorted([r for r in results if r["opt_name"] == opt_name],
                    key=lambda r: r["th"])
        if not dr:
            continue

        fname = f"{metrics_dir}/metrics_{opt_name}.txt"
        with open(fname, "w") as f:
            f.write("# threshold\tRMSE_V\tR2\tn_nonzero\n")
            for r in dr:
                f.write(f"{r['th']:.6e}\t{r['rmse_v']:.6e}\t"
                        f"{r['r2']:.6e}\t{r['n_nz']}\n")
        print(f"  Saved: {fname}")

    # Arquivo único com todos os optimizers juntos (coluna extra "optimizer")
    fname_all = f"{metrics_dir}/metrics_all_optimizers.txt"
    with open(fname_all, "w") as f:
        f.write("# optimizer\tthreshold\tRMSE_V\tR2\tn_nonzero\n")
        for r in sorted(results, key=lambda r: (r["opt_name"], r["th"])):
            f.write(f"{r['opt_name']}\t{r['th']:.6e}\t{r['rmse_v']:.6e}\t"
                    f"{r['r2']:.6e}\t{r['n_nz']}\n")
    print(f"  Saved: {fname_all}")
