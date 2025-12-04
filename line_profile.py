from input_file import InputFile
import al_fitting_functions
import numpy as np
import matplotlib.pyplot as plt
from al_fitting_functions import pseudo_voigt


class LineProfile():
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

    def __init__(self, x, y_raw):
        self.x = x
        self.y_raw = y_raw
        self.baseline = []
        self.y_corr = []
        self.fit_result = []

    def correct_y(self):
        # Baseline
        if self.baseline_method == 'asls':
            baseline = al_fitting_functions.baseline_als(self.y_raw, lam=self.asls_lambda, p=self.asls_p, niter=self.asls_niter)
        elif self.baseline_method == 'poly':
            baseline = np.polyval(np.polyfit(self.x, self.y_raw, 3), self.x)
        else:
            baseline = np.zeros_like(self.y_raw)

        # Baseline-corrected data
        y_corr = self.y_raw - baseline
        self.baseline = baseline
        self.y_corr = y_corr
    
    def fit(self, dedug:bool=False):
        if self.y_corr == []:
            self.correct_y()

        # Fit
        if dedug:
            print("Running fit (this can take some time)...")
        fit_result = al_fitting_functions.fit_spectrum(self.x, self.y_corr, self.initial_centers,
                                center_tolerance=self.center_tolerance,
                                initial_fwhm=self.initial_fwhm,
                                max_nfev=self.max_nfev,
                                verbose=False)
        if dedug:
            print("Fit result:")
            print(fit_result)
            print(fit_result['df_table'])
        self.fit_result = fit_result

    def plot_fit_and_components(self, ax, peaks_nums_2_plot:list[int]|int|None = None):
        out_prefix="fit_result"
        p_opt = self.fit_result['params']
        n_peaks = int(len(p_opt) / 4)
        # plt.figure(figsize=(10,6))


        if peaks_nums_2_plot is not None:
            if peaks_nums_2_plot == 0:
                peaks_nums_2_plot = range(n_peaks)
                ax.plot(self.x, self.y_corr, label='baseline-subtracted data', lw=1)
                ax.plot(self.x, self.fit_result['y_fit'], label='total fit', lw=1.5)
            else:
                peaks_nums_2_plot = [x - 1 for x in peaks_nums_2_plot]

            for i in peaks_nums_2_plot:
                A = p_opt[4*i + 0]; x0 = p_opt[4*i + 1]; fwhm = p_opt[4*i + 2]; eta = p_opt[4*i + 3]
                comp = pseudo_voigt(self.x, A, x0, fwhm, eta)
                ax.plot(self.x, comp, '--', lw=1, alpha=0.7, label=f'comp {i+1}' if i<6 else None)
        else:
            ax.plot(self.x, self.fit_result['y_fit'], label='total fit', lw=1.5)

        ax.legend(); plt.xlabel('x'); plt.ylabel('y'); plt.title(out_prefix + " - fit")
        # ax.tight_layout()
        # fit_png = out_prefix + ".png"
        # plt.savefig(fit_png, dpi=200); plt.close()
