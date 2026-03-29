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
def plot_critical_summary_table(df, title="Critical Analysis Summary"):
    import matplotlib.pyplot as plt

    if df is None or df.empty:
        print("[WARN] plot_critical_summary_table: boş dataframe.")
        return

    df_fmt = df.copy()

    for col in df_fmt.columns:
        if col == "System":
            continue
        df_fmt[col] = df_fmt[col].map(lambda x: f"{float(x):.2f}")

    fig, ax = plt.subplots(figsize=(10, 1.8 + 0.35 * len(df_fmt)))
    ax.axis("off")

    table = ax.table(
        cellText=df_fmt.values,
        colLabels=df_fmt.columns,
        loc="center",
        cellLoc="center"
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#D9E1F2")

    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
    plt.tight_layout()
    plt.show()
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
import numpy as np
import matplotlib.pyplot as plt

def plot_mc_with_ci(t, R, R_low, R_high):
    t = np.asarray(t, dtype=float)
    R = np.asarray(R, dtype=float)
    R_low = np.asarray(R_low, dtype=float)
    R_high = np.asarray(R_high, dtype=float)

    ci_width = R_high - R_low
    print("Max CI width:", np.max(ci_width))
    print("Mean CI width:", np.mean(ci_width))
    print("CI width at mid:", ci_width[len(ci_width)//2])

    plt.figure(figsize=(8, 5))

    plt.fill_between(
        t,
        R_low,
        R_high,
        color="orange",
        alpha=0.8,
        label="95% CI",
        zorder=1
    )

    plt.plot(
        t,
        R,
        color="blue",
        linewidth=2.5,
        label="Monte Carlo R(t)",
        zorder=3
    )

    plt.plot(
        t,
        R_low,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label="CI Lower",
        zorder=2
    )

    plt.plot(
        t,
        R_high,
        color="green",
        linestyle="--",
        linewidth=1.5,
        label="CI Upper",
        zorder=2
    )

    plt.xlabel("Time")
    plt.ylabel("Reliability R(t)")
    plt.title("Monte Carlo System Reliability with 95% CI")
    plt.ylim(
        max(0.0, np.min(R_low) - 0.02),
        min(1.02, np.max(R_high) + 0.02)
    )
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()



def plot_hazard_rate(t, R, title="Estimated Hazard Rate"):
    import numpy as np
    import matplotlib.pyplot as plt

    t = np.asarray(t, dtype=float)
    R = np.asarray(R, dtype=float)

    mask = np.isfinite(t) & np.isfinite(R)
    t = t[mask]
    R = R[mask]

    if t.size < 3:
        print("[WARN] plot_hazard_rate: yeterli veri yok.")
        return

    R_safe = np.clip(R, 1e-10, 1.0)

    dR_dt = np.gradient(R_safe, t)
    hazard = -dR_dt / R_safe

    hazard = np.nan_to_num(hazard, nan=0.0, posinf=0.0, neginf=0.0)
    hazard = np.maximum(hazard, 0.0)

    plt.figure(figsize=(7, 4))
    plt.plot(t, hazard, linewidth=2)
    plt.xlabel("Time")
    plt.ylabel("Hazard rate h(t)")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.show()
def plot_analytic_vs_mc(t_analytic, R_analytic, t_mc, R_mc, R_low=None, R_high=None):
    import numpy as np
    import matplotlib.pyplot as plt

    t_analytic = np.asarray(t_analytic, dtype=float)
    R_analytic = np.asarray(R_analytic, dtype=float)
    t_mc = np.asarray(t_mc, dtype=float)
    R_mc = np.asarray(R_mc, dtype=float)

    plt.figure(figsize=(8, 5))

    # Monte Carlo CI bandı
    if R_low is not None and R_high is not None:
        R_low = np.asarray(R_low, dtype=float)
        R_high = np.asarray(R_high, dtype=float)

        plt.fill_between(
            t_mc,
            R_low,
            R_high,
            color="orange",
            alpha=0.35,
            label="Monte Carlo 95% CI",
            zorder=1
        )

    # Analitik eğri
    plt.plot(
        t_analytic,
        R_analytic,
        color="navy",
        linewidth=2.8,
        label="Analytical / Dynamic R(t)",
        zorder=3
    )

    # Monte Carlo eğrisi
    plt.plot(
        t_mc,
        R_mc,
        color="red",
        linestyle="--",
        linewidth=2.2,
        label="Monte Carlo R(t)",
        zorder=4
    )

    plt.xlabel("Time")
    plt.ylabel("Reliability R(t)")
    plt.title("Analytical vs Monte Carlo Validation")
    plt.ylim(0, 1.05)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.show()

def build_validation_table(
    t_analytic, R_analytic,
    t_mc, R_mc,
    mttf_analytic=None,
    mttf_mc=None,
    runtime_mc=None
):
    import numpy as np
    import pandas as pd

    t_analytic = np.asarray(t_analytic, dtype=float)
    R_analytic = np.asarray(R_analytic, dtype=float)
    t_mc = np.asarray(t_mc, dtype=float)
    R_mc = np.asarray(R_mc, dtype=float)

    # MC gridinde analitik eğriyi interpolate et
    R_analytic_interp = np.interp(t_mc, t_analytic, R_analytic)

    abs_err = np.abs(R_analytic_interp - R_mc)
    rmse = np.sqrt(np.mean((R_analytic_interp - R_mc) ** 2))
    max_abs_err = np.max(abs_err)

    r_end_analytic = float(R_analytic_interp[-1])
    r_end_mc = float(R_mc[-1])
    delta_r_end = abs(r_end_analytic - r_end_mc)

    row = {
        "R(t_max) Analytic": round(r_end_analytic, 6),
        "R(t_max) Monte Carlo": round(r_end_mc, 6),
        "|ΔR(t_max)|": round(delta_r_end, 6),
        "RMSE": round(float(rmse), 6),
        "Max Abs Error": round(float(max_abs_err), 6),
    }

    if mttf_analytic is not None:
        row["MTTF Analytic"] = round(float(mttf_analytic), 4)

    if mttf_mc is not None:
        row["MTTF Monte Carlo"] = round(float(mttf_mc), 4)

    if mttf_analytic is not None and mttf_mc is not None:
        row["|ΔMTTF|"] = round(abs(float(mttf_analytic) - float(mttf_mc)), 4)

    if runtime_mc is not None:
        row["MC Runtime (s)"] = round(float(runtime_mc), 4)

    return pd.DataFrame([row])

def plot_validation_table(df, title="Analytical vs Monte Carlo Validation Summary"):
    import matplotlib.pyplot as plt

    if df is None or df.empty:
        print("[WARN] plot_validation_table: boş dataframe.")
        return

    # sayıları formatla
    df_fmt = df.copy()

    for col in df_fmt.columns:
        if "Runtime" in col:
            df_fmt[col] = df_fmt[col].map(lambda x: f"{x:.3f}")
        elif "MTTF" in col:
            df_fmt[col] = df_fmt[col].map(lambda x: f"{x:.2f}")
        else:
            df_fmt[col] = df_fmt[col].map(lambda x: f"{x:.5f}")

    fig, ax = plt.subplots(figsize=(14, 2.5))
    ax.axis("off")

    table = ax.table(
        cellText=df_fmt.values,
        colLabels=df_fmt.columns,
        loc="center",
        cellLoc="center"
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.3, 1.6)

    # header styling
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#E6E6E6")

    ax.set_title(title, fontsize=13, pad=14)

    plt.tight_layout()
    plt.show()


def plot_top_k_critical_paths(path_rts, component_paths, t_vals, k=5, title="Top-K Critical Paths at t_max"):
    import numpy as np
    import matplotlib.pyplot as plt

    if path_rts is None or len(path_rts) == 0:
        print("[WARN] plot_top_k_critical_paths: path_rts boş.")
        return

    if component_paths is None or len(component_paths) == 0:
        print("[WARN] plot_top_k_critical_paths: component_paths boş.")
        return

    # Her path için t_max noktasındaki güvenilirlik
    final_vals = []
    labels = []

    for i, r_curve in enumerate(path_rts):
        if r_curve is None or len(r_curve) == 0:
            continue

        final_r = float(r_curve[-1])
        final_vals.append(final_r)

        path_label = " - ".join(component_paths[i]) if i < len(component_paths) else f"Path {i+1}"
        labels.append(path_label)

    if not final_vals:
        print("[WARN] plot_top_k_critical_paths: çizilecek veri yok.")
        return

    # En düşük R(t_max) = en kritik yol
    order = np.argsort(final_vals)[:k]

    sorted_vals = [final_vals[i] for i in order]
    sorted_labels = [labels[i] for i in order]

    plt.figure(figsize=(10, 5))
    bars = plt.barh(range(len(sorted_vals)), sorted_vals)

    plt.yticks(range(len(sorted_vals)), sorted_labels)
    plt.xlabel("Path Reliability at t_max")
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.grid(axis="x", linestyle="--", alpha=0.4)

    for i, v in enumerate(sorted_vals):
        plt.text(v + 0.005, i, f"{v:.4f}", va="center")

    plt.tight_layout()
    plt.show()



def plot_path_contributions(path_contrib, component_paths, top_k=8, title="Monte Carlo Path Contribution"):
    import matplotlib.pyplot as plt

    if not path_contrib:
        print("[WARN] plot_path_contributions: veri yok.")
        return

    ordered = sorted(path_contrib.items(), key=lambda x: x[1], reverse=True)[:top_k]

    labels = []
    values = []

    for idx, frac in ordered:
        if idx < len(component_paths):
            path_label = " - ".join(sorted(component_paths[idx]))
        else:
            path_label = f"Path {idx+1}"

        labels.append(f"Path {idx+1}: {path_label}")
        values.append(frac * 100.0)

    plt.figure(figsize=(11, 0.6 * len(labels) + 2.2))
    bars = plt.barh(range(len(values)), values)

    plt.yticks(range(len(labels)), labels)
    plt.xlabel("Contribution to System Failure (%)")
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.grid(axis="x", linestyle="--", alpha=0.4)

    for i, v in enumerate(values):
        plt.text(v + 0.3, i, f"{v:.2f}%", va="center")

    plt.tight_layout()
    plt.show()



def plot_mc_component_importance(importance_dict, title="Monte Carlo Component Importance (ΔMTTF)"):
    import matplotlib.pyplot as plt

    if not importance_dict:
        print("[WARN] plot_mc_component_importance: veri yok.")
        return

    ordered = sorted(importance_dict.items(), key=lambda x: x[1], reverse=True)

    names = [k for k, _ in ordered]
    values = [v for _, v in ordered]

    plt.figure(figsize=(9, max(4, 0.55 * len(names) + 1.5)))
    bars = plt.barh(range(len(values)), values)

    plt.yticks(range(len(names)), names)
    plt.xlabel("ΔMTTF")
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.axvline(0, color="black", linewidth=1)
    plt.grid(axis="x", linestyle="--", alpha=0.4)

    for i, v in enumerate(values):
        x_text = v + 0.5 if v >= 0 else v - 0.5
        ha = "left" if v >= 0 else "right"
        plt.text(x_text, i, f"{v:.2f}", va="center", ha=ha)

    plt.tight_layout()
    plt.show()