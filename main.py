# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from tqdm import tqdm
import os

from input_file import InputFile
import al_IO
import al_fitting_functions


# ------------------ Load Input File -----------------------
input_file = InputFile
initial_centers = input_file.initial_centers

# ------------------ fit run options (tweak if needed) -----------------------
baseline_method = input_file.baselineOptions.baseline_method
asls_lambda = input_file.baselineOptions.asls_lambda
asls_p = input_file.baselineOptions.asls_p
asls_niter = input_file.baselineOptions.asls_niter
poly_degree = input_file.baselineOptions.poly_degree
center_tolerance = input_file.baselineOptions.center_tolerance
initial_fwhm = input_file.baselineOptions.initial_fwhm
max_nfev = input_file.baselineOptions.max_nfev

# ------------------ upload & run -------------------------------------------
filePath = r"xxxxxxx"
print("Uploaded file:", filePath)

x, y_raw = al_IO.load_xy_from_file(filePath)
print(f"Loaded {x.size} points. x range: {x.min():.3g} .. {x.max():.3g}")

# Baseline
if baseline_method == 'asls':
    baseline = al_fitting_functions.baseline_als(y_raw, lam=asls_lambda, p=asls_p, niter=asls_niter)
elif baseline_method == 'poly':
    baseline = np.polyval(np.polyfit(x, y_raw, 3), x)
else:
    baseline = np.zeros_like(y_raw)

# Baseline-corrected data
y_corr = y_raw - baseline

# Fit
print("Running fit (this can take some time)...")
fit_result = al_fitting_functions.fit_spectrum(x, y_corr, initial_centers,
                          center_tolerance=center_tolerance,
                          initial_fwhm=initial_fwhm,
                          max_nfev=max_nfev,
                          verbose=False)

# Save peak table
base = os.path.splitext(os.path.basename(filePath))[0]
out_prefix = base + "_pvfit"
df_table = fit_result['df_table']
csv_peaks = out_prefix + "_peaks.csv"
df_table.to_csv(csv_peaks, index=False)

# Build components matrix and save CSV for Excel
p_opt = fit_result['params']
n_peaks = int(len(p_opt) / 4)
components = np.zeros((x.size, n_peaks))
for i in range(n_peaks):
    A = p_opt[4*i + 0]; x0 = p_opt[4*i + 1]; fwhm = p_opt[4*i + 2]; eta = p_opt[4*i + 3]
    components[:, i] = al_fitting_functions.pseudo_voigt(x, A, x0, fwhm, eta)
y_fit = fit_result['y_fit']

# Create DataFrame with columns
out_df = pd.DataFrame({'x': x, 'y_raw': y_raw, 'baseline': baseline, 'y_corr': y_corr, 'y_fit': y_fit})
for i in range(n_peaks):
    out_df[f'comp_{i+1}'] = components[:, i]

csv_components = out_prefix + "_components.csv"
out_df.to_csv(csv_components, index=False)
print("Saved component CSV:", csv_components)

# Save plots
baseline_png, fit_png, rpng = al_IO.plot_and_save(x, y_raw, baseline, y_corr, fit_result, out_prefix=out_prefix)

print("Fit finished. Optimizer message:", fit_result.get('message', ''))
print("Saved files:")
print(" - Peak table:", csv_peaks)
print(" - Components CSV:", csv_components)
print(" - Baseline PNG:", baseline_png)
print(" - Fit PNG:", fit_png)
print(" - Residual PNG:", rpng)

# # Provide download links
# from google.colab import files as colab_files
# colab_files.download(csv_peaks)
# colab_files.download(csv_components)
# colab_files.download(baseline_png)
# colab_files.download(fit_png)
# colab_files.download(rpng)