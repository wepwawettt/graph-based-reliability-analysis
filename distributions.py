import numpy as np
import sympy
    # LogNormal hatasını çözmek için erf, erfc, sqrt eklendi
from sympy import erf, erfc, sqrt 
from sympy import symbols, expand, exp, log, Symbol
from scipy.stats import lognorm, gamma

DISTRIBUTIONS = {
    "Exponential": {
        "params": [
            {"key": "lambda", "label": "λ (Failure rate)", "type": float}
        ],
        "R_sym": lambda t, p: exp(-p["lambda"] * t),
        "sample": lambda p: np.random.exponential(1 / p["lambda"])
    },

    "Weibull": {
        "params": [
            {"key": "beta", "label": "β (Shape)", "type": float},
            {"key": "eta",  "label": "η (Scale)", "type": float}
        ],
        "R_sym": lambda t, p: exp(-(t / p["eta"]) ** p["beta"]),
        "sample": lambda p: np.random.weibull(p["beta"]) * p["eta"]
    },

    "Log-Normal": {
        "params": [
            {"key": "mu",    "label": "μ (Mean log)", "type": float},
            {"key": "sigma", "label": "σ (Std log)",  "type": float}
        ],
        "R_sym": None,  # numerik hesaplanacak
        "R_num": lambda t, p: 1 - lognorm.cdf(t, s=p["sigma"], scale=np.exp(p["mu"])),
        "sample": lambda p: np.random.lognormal(p["mu"], p["sigma"])
    },

    "Gamma": {
        "params": [
            {"key": "alpha", "label": "α (Shape)", "type": float},
            {"key": "theta", "label": "θ (Scale)", "type": float}
        ],
        "R_sym": None,  # gammaincc ile
        "R_num": lambda t, p: 1 - gamma.cdf(t, a=p["alpha"], scale=p["theta"]),
        "sample": lambda p: np.random.gamma(p["alpha"], p["theta"])
    },

    "Log-Logistic": {
        "params": [
            {"key": "alpha", "label": "α (Scale)", "type": float},
            {"key": "beta",  "label": "β (Shape)", "type": float}
        ],
        "R_sym": lambda t, p: 1 / (1 + (t / p["alpha"]) ** p["beta"]),
        "sample": lambda p: (
            p["alpha"] *
            (np.random.rand() / (1 - np.random.rand())) ** (1 / p["beta"])
        )
    },
    "Rayleigh": {
        "params": [
            {"key": "sigma", "label": "σ (Scale)", "type": float}
        ],
        "R_sym": lambda t, p: exp(-(t ** 2) / (2 * p["sigma"] ** 2)),
        "sample": lambda p: p["sigma"] * np.sqrt(-2 * np.log(np.random.rand()))
    },
    "Gompertz": {
        "params": [
            {"key": "b",   "label": "b (Shape)",  "type": float},
            {"key": "eta", "label": "η (Scale)",  "type": float}
        ],
        "R_sym": lambda t, p: exp(-p["b"] * (exp(t / p["eta"]) - 1)),
        "sample": lambda p: (
            p["eta"] * np.log(1 - np.log(np.random.rand()) / p["b"])
        )
    }
}
