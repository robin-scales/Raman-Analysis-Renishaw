import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import least_squares
import pandas as pd

# ------------------ pseudo-Voigt & model definitions ------------------------
ln2 = np.log(2.0)

def pseudo_voigt(x, A, x0, fwhm, eta):
    """Pseudo-Voigt: A = peak height, x0 center, fwhm, eta mixing fraction."""
    fwhm = np.maximum(fwhm, 1e-8)
    gauss = A * np.exp(-4*ln2 * ((x - x0)**2) / (fwhm**2))
    gamma = fwhm / 2.0
    lorentz = A * (gamma**2 / ((x - x0)**2 + gamma**2))
    return eta * lorentz + (1 - eta) * gauss

def model_sum(x, params, n_peaks):
    y = np.zeros_like(x, dtype=float)
    for i in range(n_peaks):
        A = params[4*i + 0]; x0 = params[4*i + 1]; fwhm = params[4*i + 2]; eta = params[4*i + 3]
        y += pseudo_voigt(x, A, x0, fwhm, eta)
    return y

def residuals_flat(params, x, y, n_peaks):
    return (model_sum(x, params, n_peaks) - y)

# ------------------ AsLS baseline ------------------------------------------
def baseline_als(y, lam=1e6, p=0.01, niter=10):
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, 1, 2], shape=(L-2, L))
    DTD = D.T.dot(D)
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * DTD
        z = spsolve(Z, w * y)
        residual = y - z
        w = p * (residual > 0) + (1 - p) * (residual <= 0)
    return z

# ------------------ initial params & bounds builder ------------------------
def make_initial_params(x, y, initial_centers,
                        initial_height_scale=0.2,
                        initial_fwhm=10.0,
                        initial_eta=0.5,
                        center_tolerance=8.0,
                        height_bounds=(0, np.inf),
                        fwhm_bounds=(0.5, 1000.0)):
    n_peaks = len(initial_centers)
    p0 = []; lb = []; ub = []
    ymax = np.max(y) if y.size>0 else 1.0
    # Option: set amplitude guess scaled by local density of y around center.
    for cen in initial_centers:
        A0 = ymax * initial_height_scale
        p0.append(A0); lb.append(height_bounds[0]); ub.append(height_bounds[1])
        p0.append(cen); lb.append(cen - center_tolerance); ub.append(cen + center_tolerance)
        p0.append(initial_fwhm); lb.append(fwhm_bounds[0]); ub.append(fwhm_bounds[1])
        p0.append(initial_eta); lb.append(0.0); ub.append(1.0)
    return np.array(p0, dtype=float), (np.array(lb, dtype=float), np.array(ub, dtype=float))

# ------------------ fitting routine ----------------------------------------
def fit_spectrum(x, y, initial_centers,
                 center_tolerance=8.0,
                 initial_fwhm=10.0,
                 max_nfev=5000,
                 verbose=False):
    n_peaks = len(initial_centers)
    p0, bounds = make_initial_params(x, y, initial_centers,
                                     initial_fwhm=initial_fwhm,
                                     center_tolerance=center_tolerance)
    lb, ub = bounds
    res = least_squares(residuals_flat, p0, bounds=(lb, ub),
                        args=(x, y, n_peaks),
                        max_nfev=max_nfev, verbose=2 if verbose else 0)
    p_opt = res.x
    # prepare per-peak table
    rows = []; total_area = 0.0
    for i in range(n_peaks):
        A = p_opt[4*i + 0]; x0 = p_opt[4*i + 1]; fwhm = p_opt[4*i + 2]; eta = p_opt[4*i + 3]
        window = 6.0
        xmin = x0 - window * fwhm; xmax = x0 + window * fwhm
        xi = x[(x >= xmin) & (x <= xmax)]
        if xi.size < 3:
            xi = x
        yi = pseudo_voigt(xi, A, x0, fwhm, eta)
        area = np.trapezoid(yi, xi)
        total_area += area
        rows.append({
            'Peak Index': i+1,
            'Peak Type': 'PsdVoigt1',
            'Area Intg': area,
            'FWHM': fwhm,
            'Max Height': float(np.max(yi)) if yi.size>0 else float(A),
            'Center Grvty': x0,
            'eta': eta
        })
    for r in rows:
        r['Area Intg%'] = 100.0 * (r['Area Intg'] / total_area) if total_area != 0 else 0.0
    df = pd.DataFrame(rows)
    y_fit = model_sum(x, p_opt, n_peaks)
    resids = y - y_fit
    return {'params': p_opt, 'success': res.success, 'message': res.message,
            'df_table': df, 'y_fit': y_fit, 'residuals': resids}