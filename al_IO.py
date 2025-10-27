import matplotlib.pyplot as plt
from al_fitting_functions import pseudo_voigt
import numpy as np
import pandas as pd

# ------------------ plotting + save functions ------------------------------
def plot_and_save(x, y_raw, baseline, y_corr, fit_result, out_prefix="fit_result"):
    # Baseline diagnostic
    plt.figure(figsize=(10,4))
    plt.plot(x, y_raw, label='raw', lw=1)
    plt.plot(x, baseline, label='baseline', lw=1)
    plt.plot(x, y_corr, label='baseline-subtracted', lw=1)
    plt.legend(); plt.xlabel('x'); plt.ylabel('y'); plt.title(out_prefix + " - baseline")
    plt.tight_layout()
    baseline_png = out_prefix + "_baseline.png"
    plt.savefig(baseline_png, dpi=200); plt.close()

    # Fit + components
    p_opt = fit_result['params']
    n_peaks = int(len(p_opt) / 4)
    plt.figure(figsize=(10,6))
    plt.plot(x, y_corr, label='baseline-subtracted data', lw=1)
    plt.plot(x, fit_result['y_fit'], label='total fit', lw=1.5)
    for i in range(n_peaks):
        A = p_opt[4*i + 0]; x0 = p_opt[4*i + 1]; fwhm = p_opt[4*i + 2]; eta = p_opt[4*i + 3]
        comp = pseudo_voigt(x, A, x0, fwhm, eta)
        plt.plot(x, comp, '--', lw=1, alpha=0.7, label=f'comp {i+1}' if i<6 else None)
    plt.legend(); plt.xlabel('x'); plt.ylabel('y'); plt.title(out_prefix + " - fit")
    plt.tight_layout()
    fit_png = out_prefix + ".png"
    plt.savefig(fit_png, dpi=200); plt.close()

    # Residuals
    plt.figure(figsize=(10,2.5))
    plt.plot(x, fit_result['residuals'], label='residuals')
    plt.axhline(0, color='k', lw=0.5)
    plt.tight_layout()
    rname = out_prefix + "_residuals.png"
    plt.savefig(rname, dpi=200); plt.close()

    return baseline_png, fit_png, rname

# ------------------ loader for .txt/.csv with auto-delimiter ----------------
def load_xy_from_file(path):
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if s == "" or s.startswith('#'):
                continue
            first_line = s
            break
        else:
            raise ValueError("File appears empty or contains only comments.")
    if ',' in first_line:
        sep = ','
    elif '\t' in first_line:
        sep = '\t'
    else:
        sep = r'\s+'
    df = pd.read_csv(path, sep=sep, header=None, comment='#', engine='python')
    if df.shape[1] < 2:
        raise ValueError("File must contain at least two columns (x y).")
    x = df.iloc[:,0].astype(float).values
    y = df.iloc[:,1].astype(float).values
    order = np.argsort(x)
    return x[order], y[order]