# Graph-Based System Reliability Analysis Tool

This repository contains my undergraduate graduation project developed at **Ankara University – Computer Engineering Department**.  
The project provides a **graph-based analytical and simulation-driven reliability analysis tool** for complex multi-component systems.

The system supports **static reliability**, **time-dependent reliability R(t)**, and **Monte Carlo simulation**, all integrated into an interactive **PyQt6-based GUI**.

---

## 📌 Project Overview

Modern engineering systems consist of interconnected components whose failures are often interdependent.  
This project models such systems as **graphs** and evaluates system reliability using:

- Analytical reliability theory
- Inclusion–exclusion principle
- Monte Carlo simulation
- Path-based and component-based criticality analysis

The tool enables both **quantitative reliability evaluation** and **qualitative insight** into critical components and degradation behavior.

---

## 🚀 Key Features

- Interactive graph-based system modeling
- Automatic extraction of **minimal path sets**
- Static system reliability analysis
- Time-dependent reliability analysis \( R(t) \)
- Supported lifetime distributions:
  - Exponential
  - Weibull
  - Log-Normal
  - Gamma
  - Log-Logistic
  - Rayleigh
  - Gompertz
- Monte Carlo simulation using analytical path sets
- Common Cause Failure (CCF) modeling using β-factor
- System lifetime histogram
- Path reliability curves
- Critical interval detection (t₉₀ – t₁₀)
- Component Criticality Index (CCI)
- Sensitivity analysis (Tornado chart)
- Multi-model critical comparison
- Fully integrated PyQt6 graphical interface

---

## 🧠 Methodology

1. The system is modeled as a graph with components and junctions.
2. All minimal paths between **Start** and **End** nodes are extracted.
3. System reliability is derived analytically using the **inclusion–exclusion principle**.
4. Time-dependent reliability \( R(t) \) is computed from component lifetime distributions.
5. Monte Carlo simulation generates system lifetime samples using the same path sets.
6. Analytical and simulation-based results are compared and validated.
7. Criticality and robustness analyses are performed at system, path, and component levels.

> Monte Carlo simulation is intentionally performed **after analytical path extraction** to ensure methodological consistency.

---

## 📊 Monte Carlo Simulation

Monte Carlo simulation produces:
- System lifetime samples
- Empirical reliability function:

\[
R_{MC}(t) = P(T_{system} > t)
\]

These results are used to validate analytical reliability curves and analyze uncertainty.

---

## 🖥️ Project Structure

```text
├── main.py                # GUI and application logic
├── distributions.py       # Lifetime distributions and sampling functions
├── monte_carlo.py         # Monte Carlo simulation engine
├── critical_analysis.py   # Criticality, robustness and sensitivity analysis
├── README.md
⚙️ Installation
Requirements
Python 3.10 or higher

NumPy

SciPy

SymPy

Matplotlib

PyQt6

Install dependencies:

bash
Kodu kopyala
pip install numpy scipy sympy matplotlib pyqt6
▶️ Usage
Run the application:

bash
Kodu kopyala
python main.py
Typical workflow:

Create or load a system model

Select analysis mode:

Static Reliability

Dynamic Reliability R(t)

Monte Carlo Simulation

Define component parameters and connections

Run analysis

Visualize results and criticality metrics

📈 Example Outputs
System reliability curve 
𝑅
(
𝑡
)
R(t)

Monte Carlo reliability curve

System lifetime histogram

Path reliability evolution

Critical interval plots

Component Criticality Index (CCI)

Sensitivity tornado charts

🔬 Academic Context
This project was developed as an undergraduate graduation project and follows standard reliability engineering methodologies.

It is intended for:

Academic use

Educational demonstrations

Research-oriented extensions

🔮 Future Work
Bayesian reliability modeling

Graph Neural Networks (GNNs) for reliability prediction

Large-scale system optimization

Confidence interval estimation and uncertainty quantification

👩‍💻 Author
Selin Ayhan
Computer Engineering
Ankara University

📄 License
This project is provided for academic and educational purposes only.
