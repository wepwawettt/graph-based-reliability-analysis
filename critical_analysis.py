import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker

def find_crossing_time(t, R, level):
    t = np.asarray(t, dtype=float)
    R = np.asarray(R, dtype=float)

    mask = np.isfinite(t) & np.isfinite(R)
    t = t[mask]
    R = R[mask]

    if np.all(R > level):
        return None

    k = np.argmax(R <= level)

    if k == 0:
        return t[0]

    t1, t2 = t[k-1], t[k]
    r1, r2 = R[k-1], R[k]

    return t1 + (level - r1) * (t2 - t1) / (r2 - r1)


def plot_critical_intervals(results_dict):
    intervals = {}

    for name, data in results_dict.items():
        t = data["t"]
        R = data["R"]

        t90 = find_crossing_time(t, R, 0.9)
        t10 = find_crossing_time(t, R, 0.1)

        if t90 is None or t10 is None:
            continue

        intervals[name] = (t90, t10)

    fig, ax = plt.subplots(figsize=(8, 4))

    for i, (name, (t90, t10)) in enumerate(intervals.items()):
        ax.hlines(i, t90, t10, linewidth=6)
        t_max_local = np.max(results_dict[name]["t"])
        ax.vlines([t90, t10], i - 0.15, i + 0.15, linestyles="dashed")
        ax.text(t90, i + 0.2, f"{t90:.1f}", ha="center", fontsize=9)
        ax.text(t10, i + 0.2, f"{t10:.1f}", ha="center", fontsize=9)

    ax.set_yticks(range(len(intervals)))
    ax.set_yticklabels(intervals.keys())
    ax.set_xlabel("Time")
    ax.set_title("Critical Reliability Window (0.9 – 0.1)")
    plt.tight_layout()
    plt.show()

    return intervals

def plot_critical_slopes(intervals):
    slopes = {}

    for name, (t90, t10) in intervals.items():
        if t10 > t90:
            slopes[name] = (0.9 - 0.1) / (t10 - t90)
        else:
            slopes[name] = np.nan

    names = list(slopes.keys())
    values = list(slopes.values())

    plt.figure(figsize=(6, 4))
    bars = plt.bar(names, values)

    plt.ylabel("Average Degradation Rate")
    plt.title("Critical Region Degradation Rate")

    # 👉 SAYILARI BAR ÜSTÜNE YAZ
    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{value:.6f}",   # 👈 6 basamak
            ha="center",
            va="bottom",
            fontsize=9
        )

    plt.tight_layout()
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.6f'))
    plt.show()

def critical_summary_table(intervals):
    rows = []

    for name, (t90, t10) in intervals.items():
        rows.append({
            "System": name,
            "t_90": round(t90, 2),
            "t_10": round(t10, 2),
            "Critical Duration": round(t10 - t90, 2)
        })

    return pd.DataFrame(rows)

import matplotlib.pyplot as plt
import numpy as np

def plot_component_criticality(results_dict, paths_dict):
    """
    results_dict:
      {model: {"t":[], "R":[]} }

    paths_dict:
      {model: [set(a1,a2), set(a3,a4)...]}
    """
    model_name = list(results_dict.keys())[0]

    component_scores = {}

    for model, paths in paths_dict.items():
        t = np.array(results_dict[model]["t"])
        R = np.array(results_dict[model]["R"])

        # t90, t10
        t90 = t[np.argmax(R <= 0.9)]
        t10 = t[np.argmax(R <= 0.1)]
        delta_t = max(t10 - t90, 1e-6)

        for path in paths:
            for comp in path:
                component_scores.setdefault(comp, 0)
                component_scores[comp] += 1 / (len(path) * delta_t)

    names = list(component_scores.keys())
    values = list(component_scores.values())

    plt.figure(figsize=(7, 4))
    bars = plt.bar(names, values)

    plt.title(f"Component Criticality Index (CCI) – {model_name}")

    plt.ylabel("Criticality Score")
    plt.xlabel("Component")

    for bar, val in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height(),
            f"{val:.2e}",          # ✅ .4f yerine .2e
            ha="center", va="bottom", fontsize=9
        )


    plt.tight_layout()
    plt.show()
    import matplotlib.ticker as mticker
    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))

def plot_path_robustness(results_dict, paths_dict):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(figsize=(6, 4))

    for model, paths in paths_dict.items():
        t = np.array(results_dict[model]["t"])
        R = np.array(results_dict[model]["R"])

        t90 = t[np.argmax(R <= 0.9)]
        t10 = t[np.argmax(R <= 0.1)]
        delta_t = t10 - t90

        plt.scatter(t90, delta_t, s=80, label=model)

    plt.xlabel("t₉₀ (Start of degradation)")
    plt.ylabel("Δt = t₁₀ − t₉₀ (Critical duration)")
    plt.title("Path Robustness Comparison")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_sensitivity_tornado(mttf_base, sensitivity_results):
    """
    sensitivity_results = {
      "a1": +120,
      "a2": -80,
      ...
    }
    """

    names = list(sensitivity_results.keys())
    effects = list(sensitivity_results.values())

    plt.figure(figsize=(6, 4))
    plt.barh(names, effects)
    plt.axvline(0, color="black")

    plt.xlabel("ΔMTTF")
    plt.title("Sensitivity Analysis (Tornado Chart)")
    plt.tight_layout()
    plt.show()


import numpy as np
import matplotlib.pyplot as plt

def plot_true_path_robustness(results_dict, paths_dict, path_rt_dict):
    """
    results_dict:
        {model: {"t": [...], "R": [...]}}

    paths_dict:
        {model: [set(a1,a2), set(a3,a4), ...]}

    path_rt_dict:
        {model: [R_path1(t), R_path2(t), ...]}
    """

    markers = ["o", "s", "^", "D", "v", "P", "X"]
    plt.figure(figsize=(7, 5))

    for m_idx, model in enumerate(results_dict):
        t = np.array(results_dict[model]["t"])
        paths = paths_dict[model]
        path_rts = path_rt_dict[model]

        for i, (path, R_path) in enumerate(zip(paths, path_rts)):
            R_path = np.array(R_path)

            # t90, t10 (path bazlı)
            if np.all(R_path > 0.9) or np.all(R_path > 0.1):
                continue

            t90 = t[np.argmax(R_path <= 0.9)]
            t10 = t[np.argmax(R_path <= 0.1)]
            delta_t = t10 - t90

            label = f"{model} | Path {i+1}"

            plt.scatter(
                t90,
                delta_t,
                marker=markers[i % len(markers)],
                s=90,
                label=label
            )

    plt.xlabel("t₉₀ (Path degradation start)")
    plt.ylabel("Δt = t₁₀ − t₉₀ (Path critical duration)")
    plt.title("True Path Robustness Comparison")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.show()


def plot_path_rt_curves(results_dict, paths_dict, path_rt_dict):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(figsize=(8, 5))

    for model in results_dict:
        t = np.array(results_dict[model]["t"])
        paths = paths_dict[model]
        path_rts = path_rt_dict[model]

        for i, (path, R_path) in enumerate(zip(paths, path_rts)):
            label = f"{model} | Path {i+1}: {'-'.join(sorted(path))}"

            plt.plot(
                t,
                R_path,
                linewidth=2,
                label=label
            )

    plt.axhline(0.9, linestyle="--", color="gray", alpha=0.6)
    plt.axhline(0.1, linestyle="--", color="gray", alpha=0.6)

    plt.xlabel("Time")
    plt.ylabel("Path Reliability Rₚ(t)")
    plt.title("Path Reliability Evolution Over Time")
    plt.legend(fontsize=7)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()


# ===== EK: Survival + CDF (Histogram samples üzerinden) =====
import numpy as np
import matplotlib.pyplot as plt

def plot_survival_and_cdf(T_sys_samples, title_prefix="System Lifetime"):
    """
    T_sys_samples: 1D array-like (Monte Carlo / simülasyon sistem ömrü örnekleri)
    Survival: S(t)=P(T>t)  (bu aslında R(t) ile aynı anlamda kullanılabilir)
    CDF: F(t)=P(T<=t)
    """
    T = np.asarray(T_sys_samples, dtype=float)
    T = T[np.isfinite(T)]
    if T.size < 5:
        print("[WARN] plot_survival_and_cdf: Yeterli sample yok.")
        return

    T_sorted = np.sort(T)
    n = T_sorted.size

    # Empirical CDF: F(t_i)=i/n
    F = np.arange(1, n + 1) / n
    # Survival: S(t)=1-F(t)
    S = 1.0 - F

    # --- CDF ---
    plt.figure(figsize=(7, 4))
    plt.plot(T_sorted, F, linewidth=2)
    plt.xlabel("t")
    plt.ylabel("F(t) = P(T ≤ t)")
    plt.title(f"{title_prefix} - Empirical CDF")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()

    # --- Survival ---
    plt.figure(figsize=(7, 4))
    plt.plot(T_sorted, S, linewidth=2)
    plt.xlabel("t")
    plt.ylabel("S(t) = P(T > t)")
    plt.title(f"{title_prefix} - Survival Function (Reliability)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()


import os

def plot_system_rt_comparison(results_dict):
    plt.figure(figsize=(8, 5))

    for name, data in results_dict.items():
        t = np.asarray(data["t"])
        R = np.asarray(data["R"])

        clean_name = os.path.splitext(name)[0]   # 👈 BURASI
        plt.plot(t, R, linewidth=2, label=clean_name)

    plt.axhline(0.9, linestyle="--", color="gray", alpha=0.5)
    plt.axhline(0.1, linestyle="--", color="gray", alpha=0.5)

    plt.xlabel("Zaman (t)")
    plt.ylabel("Sistem Güvenilirliği R(t)")
    plt.title("Sistem Bazlı Güvenilirlik Eğrilerinin Karşılaştırılması")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()
