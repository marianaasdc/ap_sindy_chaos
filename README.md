# Data-Driven Discovery and Reduction of Action Potential Dynamics with SINDy

---

## Overview

<!--This repository contains the code and notebooks associated with the paper:

> *Data-driven discovery and reduction of action potential dynamics with SINDy*
> M. A. S. de Carvalho, R. W. dos Santos, and B. M. Rocha *(submitted)* 
-->

Code that uses **SINDy** (Sparse Identification of Nonlinear Dynamics) to automatically discover and/or reduce the governing equations of electrophysiological models from simulated time-series data. The method is evaluated on the Hodgkin-Huxley (HH) and FitzHugh-Nagumo (FHN) models, including PCA-based dimensionality reduction and constrained identification.

---

## Repository Structure

```
ap_sindy_chaos/
│
├── model_fhn.py                     # FitzHugh-Nagumo ODE model & solver
├── model_hh.py                      # Hodgkin-Huxley ODE model & solver
│
├── sindy_fhn.ipynb                  # SINDy identification on FHN model
├── sindy_hh_contrained.ipynb        # Constrained SINDy on HH model
├── sindy_hh_pca.ipynb               # SINDy on PCA-reduced HH data
├── sindy_hh_pca_fhnlike.ipynb       # SINDy on PCA-reduced HH (FHN-like regime)
│
├── test_model_fhn.ipynb             # FHN model simulation & validation
├── test_model_hh.ipynb              # HH model simulation & validation
├── test_sindy_hh_constrained.ipynb  # Tests for constrained SINDy
└── test_sindy_hh_pca.ipynb          # Tests for PCA + SINDy pipeline
```

---

## Methods

Simulated trajectories are generated from the HH and FHN models and used as input to SINDy via the [PySINDy](https://github.com/dynamicslab/pysindy) library. Several sparse regression optimizers are benchmarked:

- **SR3** — Sparse Relaxed Regularized Regression
- **STLSQ** — Sequentially Thresholded Least Squares
- **LASSO** — L1-regularized regression
- **SSR** — Stepwise Sparse Regression


---

## Requirements

```bash
pip install pysindy numpy scipy matplotlib scikit-learn
```

| Package | Version |
|---|---|
| Python | 3.11 |
| pysindy | 1.7.5 |
| numpy | 1.24.4 |
| scipy | 1.11.0 |
| matplotlib | 3.7.0 |
| scikit-learn | 1.6.1 |
| pandas | 1.5.3 |
| sympy | 1.13.1 |
| jupyter | 1.1.1 |


---

## Authors 

[M. A. S. de Carvalho](https://github.com/marianaasdc/), R. W. dos Santos, and [B. M. Rocha](https://github.com/rochabm)



<!--
## Citation

If you use this code, please cite:

```
M. A. S. de Carvalho, R. W. dos Santos, and B. M. Rocha,
"Data-driven discovery and reduction of action potential dynamics with SINDy", (submitted).
```
-->