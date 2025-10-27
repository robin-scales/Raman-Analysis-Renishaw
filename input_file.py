# ------------------ Input File ------------------------

# https://docs.python.org/3/tutorial/classes.html
class BaselineOptions:
    """ Baseline options (tweak if needed)
    """
    # Methods supported: 'asls' (Asymmetric Least Squares), 'poly' (polynomial fit), or None
    baseline_method = 'asls'   # default: 'asls'
    
    # AsLS params (good starting values for Raman)
    asls_lambda = 1e6   # smoothing; increase for smoother baseline (e.g., 1e7)
    asls_p = 0.01       # asymmetry (0<p<1) smaller p treats peaks as positive; typical 0.001-0.1
    asls_niter = 10     # iterations (8-20)
    
    center_tolerance = 8
    initial_fwhm = 5.0
    max_nfev = 5000

    # Polynomial baseline params (only if baseline_method == 'poly')
    poly_degree = 3

class InputFile:
    """ Baseline options (tweak if needed)
    """
    material_name = "GaN-on-Si"

    baselineOptions = BaselineOptions

    # ------------------ initial centers (user-provided) -------------------------
    # initial_centers = [
    #     521.95, 567, 734
    # ]
    initial_centers = [
        521.95, 566.273, 734
    ]

    peak_names = [
        "Si", "GaN E2 (high)", "GaN A1 (LO)"
    ]

    resouces = ["https://doi.org/10.1002/sia.1134", "https://ramanlife.com/library/gallium-nitride/", "https://ramanlife.com/library/silicon-raman/"]

    # Allow any number, use len(initial_centers) dynamically
    n_peaks_global = len(initial_centers)
    assert n_peaks_global > 0, "Provide at least one center in initial_centers."
