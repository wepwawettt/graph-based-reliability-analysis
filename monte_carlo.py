import numpy as np
import matplotlib.pyplot as plt
import copy
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# 1) Tek bileşen için lifetime örnekleme
#    Burada local rng kullanıyoruz -> seed gerçekten çalışsın
# =========================================================
def sample_one_lifetime(dist_name, params, rng):
    if dist_name == "static":
        return 1e20

    if dist_name == "Exponential":
        lam = float(params["lambda"])
        return rng.exponential(1.0 / lam)

    if dist_name == "Weibull":
        beta = float(params["beta"])
        eta = float(params["eta"])
        return rng.weibull(beta) * eta

    if dist_name == "Log-Normal":
        mu = float(params["mu"])
        sigma = float(params["sigma"])
        return rng.lognormal(mean=mu, sigma=sigma)

    if dist_name == "Gamma":
        alpha = float(params["alpha"])
        theta = float(params["theta"])
        return rng.gamma(shape=alpha, scale=theta)

    if dist_name == "Log-Logistic":
        alpha = float(params["alpha"])
        beta = float(params["beta"])
        u = np.clip(rng.random(), 1e-12, 1 - 1e-12)
        return alpha * (u / (1.0 - u)) ** (1.0 / beta)

    if dist_name == "Rayleigh":
        sigma = float(params["sigma"])
        return rng.rayleigh(scale=sigma)

    if dist_name == "Gompertz":
        b = float(params["b"])
        eta = float(params["eta"])
        u = np.clip(rng.random(), 1e-12, 1 - 1e-12)
        return eta * np.log(1.0 - np.log(1.0 - u) / b)

    raise ValueError(f"Bilinmeyen dağılım tipi: {dist_name}")


# =========================================================
# 2) Tüm bileşenler için lifetime üret
# =========================================================
def sample_component_lifetimes(components, rng, ccf=None):
    lifetimes = {}
    T_ccf = None

    if ccf:
        beta, lambdas = ccf
        if beta > 0 and lambdas:
            lambda_avg = np.mean(lambdas)
            T_ccf = rng.exponential(1.0 / (beta * lambda_avg))

    for cname, d in components.items():
        lt_ind = sample_one_lifetime(d["dist"], d.get("params", {}), rng)

        if T_ccf is not None:
            lifetimes[cname] = min(lt_ind, T_ccf)
        else:
            lifetimes[cname] = lt_ind

    return lifetimes


# =========================================================
# 3) Path fail times
# =========================================================
def compute_path_fail_times(lifetimes, component_paths):
    return [
        min(lifetimes[c] for c in path)
        for path in component_paths if path
    ]


# =========================================================
# 4) System failure time
# =========================================================
def compute_system_failure_time(path_fail_times):
    return max(path_fail_times) if path_fail_times else 1e20


# =========================================================
# 5) R(t) eğrisi
# =========================================================
def estimate_reliability_curve(T_sys, t_max, n_t):
    t_vals = np.linspace(0, t_max, n_t)
    R_mc = np.array([np.mean(T_sys > t) for t in t_vals], dtype=float)
    return t_vals, R_mc


# =========================================================
# 6) Güven aralığı
# =========================================================
def compute_reliability_ci(R_mc, N):
    se = np.sqrt(np.maximum(R_mc * (1.0 - R_mc), 0.0) / N)
    r_low = np.clip(R_mc - 1.96 * se, 0.0, 1.0)
    r_high = np.clip(R_mc + 1.96 * se, 0.0, 1.0)
    return r_low, r_high


# =========================================================
# 7) MTTF + CI
# =========================================================
def compute_mttf_stats(T_sys):
    T_sys = np.asarray(T_sys, dtype=float)
    n = len(T_sys)

    if n == 0:
        return np.nan, np.nan, np.nan

    mttf = np.mean(T_sys)

    if n == 1:
        return float(mttf), float(mttf), float(mttf)

    std = np.std(T_sys, ddof=1)
    half_width = 1.96 * std / np.sqrt(n)

    ci_low = mttf - half_width
    ci_high = mttf + half_width
    return float(mttf), float(ci_low), float(ci_high)


# =========================================================
# 8) Ana Monte Carlo
# =========================================================
def run_monte_carlo(components, component_paths, N, t_max, ccf=None, n_t=100, seed=None):
    rng = np.random.default_rng(seed)
    T_sys = np.zeros(N, dtype=float)

    # her path'in sisteme katkı sayacı
    path_contrib_counts = np.zeros(len(component_paths), dtype=float)

    for i in range(N):
        lifetimes = sample_component_lifetimes(components, rng, ccf=ccf)
        path_fail_times = compute_path_fail_times(lifetimes, component_paths)
        T_sys[i] = compute_system_failure_time(path_fail_times)

        # sistemi belirleyen yol(lar): en geç fail olan path
        if path_fail_times:
            max_ft = max(path_fail_times)
            winners = [
                j for j, ft in enumerate(path_fail_times)
                if np.isclose(ft, max_ft, rtol=1e-9, atol=1e-12)
            ]

            if winners:
                share = 1.0 / len(winners)
                for j in winners:
                    path_contrib_counts[j] += share

    t_vals, R_mc = estimate_reliability_curve(T_sys, t_max, n_t)
    R_low, R_high = compute_reliability_ci(R_mc, N)
    MTTF, CI_low, CI_high = compute_mttf_stats(T_sys)

    path_contrib = {
        idx: float(cnt / N)
        for idx, cnt in enumerate(path_contrib_counts)
        if cnt > 0
    }

    return T_sys, t_vals, R_mc, R_low, R_high, MTTF, CI_low, CI_high, path_contrib


# =========================================================
# 9) Convergence analizi için checkpoint seç
# =========================================================
def build_convergence_checkpoints(N, min_points=10, max_points=25):
    if N <= 1:
        return np.array([1], dtype=int)

    start = min(min_points, N)
    raw = np.geomspace(start, N, num=min(max_points, N)).astype(int)
    raw = np.unique(raw)

    if raw[0] != 1:
        raw = np.insert(raw, 0, 1)
    if raw[-1] != N:
        raw = np.append(raw, N)

    return raw


# =========================================================
# 10) Cumulative MTTF convergence
# =========================================================
def monte_carlo_convergence(T_sys, analytic_mttf=None, title="Monte Carlo Convergence (MTTF)"):
    T_sys = np.asarray(T_sys, dtype=float)
    T_sys = T_sys[np.isfinite(T_sys)]

    if T_sys.size < 2:
        print("[WARN] monte_carlo_convergence: yeterli sample yok.")
        return None

    n_total = len(T_sys)
    checkpoints = build_convergence_checkpoints(n_total)

    csum = np.cumsum(T_sys)
    csum_sq = np.cumsum(T_sys ** 2)

    means = []
    ci_lows = []
    ci_highs = []

    for n in checkpoints:
        mean_n = csum[n - 1] / n

        if n == 1:
            var_n = 0.0
        else:
            ex2 = csum_sq[n - 1] / n
            var_n = max(ex2 - mean_n ** 2, 0.0) * n / max(n - 1, 1)

        std_n = np.sqrt(var_n)
        half_width = 1.96 * std_n / np.sqrt(n)

        means.append(mean_n)
        ci_lows.append(mean_n - half_width)
        ci_highs.append(mean_n + half_width)

    means = np.asarray(means, dtype=float)
    ci_lows = np.asarray(ci_lows, dtype=float)
    ci_highs = np.asarray(ci_highs, dtype=float)

    plt.figure(figsize=(8, 5))
    plt.plot(checkpoints, means, linewidth=2.2, label="Monte Carlo cumulative MTTF")
    plt.fill_between(checkpoints, ci_lows, ci_highs, alpha=0.25, label="95% CI")

    if analytic_mttf is not None and np.isfinite(analytic_mttf):
        plt.axhline(
            analytic_mttf,
            linestyle="--",
            linewidth=2,
            label=f"Analytical MTTF = {analytic_mttf:.2f}"
        )

    plt.xlabel("Simulation Count (N)")
    plt.ylabel("Estimated MTTF")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()

    final_mean = means[-1]
    rel_change = None
    if len(means) >= 2 and abs(final_mean) > 1e-12:
        rel_change = abs(means[-1] - means[-2]) / abs(final_mean)

    print("===== Monte Carlo Convergence Summary =====")
    print("Final sample count :", n_total)
    print("Final MTTF estimate:", round(float(final_mean), 6))
    if rel_change is not None:
        print("Last relative change:", f"{rel_change:.6e}")

    return {
        "N_points": checkpoints,
        "mean_mttf": means,
        "ci_low": ci_lows,
        "ci_high": ci_highs,
        "analytic_mttf": analytic_mttf,
        "final_mttf": float(final_mean),
        "last_relative_change": None if rel_change is None else float(rel_change),
    }

def _improve_component_for_importance(comp_data, delta=0.10):
    d = copy.deepcopy(comp_data)

    if d["dist"] == "static":
        if "R" in d:
            d["R"] = min(0.999999, d["R"] + (1.0 - d["R"]) * delta)
        return d

    p = d.get("params", {})
    dist = d["dist"]

    if dist == "Exponential" and "lambda" in p:
        p["lambda"] = max(1e-12, p["lambda"] * (1.0 - delta))

    elif dist == "Weibull" and "eta" in p:
        p["eta"] *= (1.0 + delta)

    elif dist == "Log-Normal" and "mu" in p:
        p["mu"] += np.log(1.0 + delta)

    elif dist == "Gamma" and "theta" in p:
        p["theta"] *= (1.0 + delta)

    elif dist == "Log-Logistic" and "alpha" in p:
        p["alpha"] *= (1.0 + delta)

    elif dist == "Rayleigh" and "sigma" in p:
        p["sigma"] *= (1.0 + delta)

    elif dist == "Gompertz":
        if "eta" in p:
            p["eta"] *= (1.0 + delta)
        if "b" in p:
            p["b"] = max(1e-12, p["b"] * (1.0 - 0.5 * delta))

    return d


def monte_carlo_component_importance(
    components,
    component_paths,
    N,
    t_max,
    ccf=None,
    delta=0.10,
    seed=42
):
    _, _, _, _, _, base_mttf, _, _, _ = run_monte_carlo(
        components=components,
        component_paths=component_paths,
        N=N,
        t_max=t_max,
        ccf=ccf,
        seed=seed
    )

    importance = {}

    for cname in components:
        modified = copy.deepcopy(components)
        modified[cname] = _improve_component_for_importance(modified[cname], delta=delta)

        _, _, _, _, _, new_mttf, _, _, _ = run_monte_carlo(
            components=modified,
            component_paths=component_paths,
            N=N,
            t_max=t_max,
            ccf=ccf,
            seed=seed
        )

        importance[cname] = float(new_mttf - base_mttf)

    return float(base_mttf), importance