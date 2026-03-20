import numpy as np
from distributions import DISTRIBUTIONS


def run_monte_carlo(components, component_paths, N, t_max, ccf=None, n_t=100, seed=None):
    """
    Monte Carlo simulation with confidence interval and MTTF statistics

    Returns:
        T_sys   : system lifetime samples
        t_vals  : time grid
        R_mc    : Monte Carlo reliability estimate
        R_lower : lower 95% CI for R(t)
        R_upper : upper 95% CI for R(t)
        MTTF    : mean time to failure
        CI_low  : lower 95% CI for MTTF
        CI_high : upper 95% CI for MTTF
    """

    rng = np.random.default_rng(seed)

    T_sys = np.zeros(N, dtype=float)
    
    for i in range(N):
        lifetimes = {}
        T_ccf = None

        # ----- CCF event -----
        if ccf:
            beta, lambdas = ccf
            if beta > 0 and lambdas:
                lambda_avg = np.mean(lambdas)
                T_ccf = rng.exponential(1 / (beta * lambda_avg))

        # ----- component lifetimes -----
        for cname, d in components.items():
            if d["dist"] == "static":
                lt_ind = 1e20
            else:
                conf = DISTRIBUTIONS[d["dist"]]
                lt_ind = conf["sample"](d["params"])

            lifetimes[cname] = min(lt_ind, T_ccf) if T_ccf is not None else lt_ind

        # ----- path failures -----
        path_fail_times = [
            min(lifetimes[c] for c in p)
            for p in component_paths if p
        ]

        T_sys[i] = max(path_fail_times) if path_fail_times else 1e20

    # ----- reliability curve -----
    t_vals = np.linspace(0, t_max, n_t)

    R_mc = np.array([
        np.mean(T_sys > t) for t in t_vals
    ], dtype=float)

    # ----- confidence interval for R(t) -----
    SE = np.sqrt(R_mc * (1 - R_mc) / N)

    R_lower = np.clip(R_mc - 1.96 * SE, 0, 1)
    R_upper = np.clip(R_mc + 1.96 * SE, 0, 1)

    # ----- MTTF statistics -----
    MTTF = np.mean(T_sys)
    STD = np.std(T_sys, ddof=1)

    CI_low = MTTF - 1.96 * STD / np.sqrt(N)
    CI_high = MTTF + 1.96 * STD / np.sqrt(N)

    
    return (
        T_sys,
        t_vals,
        R_mc,
        R_lower,
        R_upper,
        MTTF,
        CI_low,
        CI_high
    )
    

def monte_carlo_convergence(T_sys, analytic_mttf=None, title="Monte Carlo Convergence Study"):
    import numpy as np
    import matplotlib.pyplot as plt

    if T_sys is None or len(T_sys) < 2:
        print("[WARN] monte_carlo_convergence: yetersiz T_sys.")
        return

    # küçük N'lerden başlayıp tam N'e kadar gidelim
    n_min = min(50, len(T_sys))
    if len(T_sys) < 50:
        n_min = 2

    N_vals = np.linspace(n_min, len(T_sys), 20, dtype=int)
    N_vals = np.unique(N_vals)

    means = []
    ci_lows = []
    ci_highs = []

    for n in N_vals:
        sample = T_sys[:n]
        mean_n = np.mean(sample)
        std_n = np.std(sample, ddof=1) if n > 1 else 0.0
        se_n = std_n / np.sqrt(n) if n > 1 else 0.0
        ci = 1.96 * se_n

        means.append(mean_n)
        ci_lows.append(mean_n - ci)
        ci_highs.append(mean_n + ci)

    plt.figure(figsize=(8, 5))
    plt.plot(N_vals, means, marker="o", label="MC Estimated MTTF")
    plt.fill_between(N_vals, ci_lows, ci_highs, alpha=0.2, label="95% CI")

    if analytic_mttf is not None:
        plt.axhline(
    y=analytic_mttf,
    linestyle="--",
    linewidth=3,
    color="red"
)
    plt.xlabel("Sample Size (N)")
    plt.ylabel("Estimated MTTF")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()