# Graph-Based System Reliability Analysis Tool

A graph-based reliability analysis framework for multi-component systems, developed as an undergraduate graduation project at **Ankara University – Computer Engineering Department**.

This project combines **analytical reliability evaluation**, **time-dependent reliability analysis**, and **Monte Carlo simulation** in an interactive **PyQt6-based desktop application**. It is designed to support both **system-level reliability estimation** and **path/component-level criticality interpretation**.

---

## Overview

Modern engineering systems often consist of interconnected components with heterogeneous failure behaviors.  
This project models such systems as graphs and evaluates reliability using:

- graph-based system representation
- minimal path extraction
- inclusion-exclusion-based analytical reliability computation
- dynamic reliability analysis \(R(t)\)
- Monte Carlo validation
- criticality and robustness analysis

The tool is intended not only to compute reliability metrics, but also to help explain **which paths and components are most critical** in the system.

---

## Key Features

- Interactive graph-based system modeling
- Automatic extraction of **minimal path sets**
- **Static reliability analysis**
- **Dynamic reliability analysis** with time-dependent \(R(t)\)
- Support for multiple component lifetime distributions:
  - Exponential
  - Weibull
  - Log-Normal
  - Gamma
  - Log-Logistic
  - Rayleigh
  - Gompertz
- **Monte Carlo simulation**
- Monte Carlo **95% confidence interval** visualization
- **Analytical vs Monte Carlo validation**
- Validation summary metrics:
  - RMSE
  - maximum absolute error
  - terminal reliability difference
  - MTTF difference
  - runtime
- **Monte Carlo convergence analysis**
- **Hazard rate analysis**
- **Critical interval analysis** \((t_{90}, t_{10})\)
- **Component Criticality Index (CCI)**
- **Sensitivity / tornado analysis**
- **Top-k critical path analysis**
- **Path contribution analysis**
- **Monte Carlo component importance**
- Optional **Common Cause Failure (CCF)** modeling
- Model **save/load** support
- Multi-model **comparison and critical analysis**
- Interactive **PyQt6 GUI**

---

## Methodology

The framework follows the steps below:

1. Model the system as a graph between **Start** and **End** nodes  
2. Extract minimal component paths  
3. Build system reliability using inclusion-exclusion logic  
4. Compute analytical or dynamic reliability curves  
5. Validate results with Monte Carlo simulation  
6. Interpret vulnerability through path/component criticality analyses  

---

## Monte Carlo Reliability

The empirical system reliability is estimated as:

\[
R_{MC}(t) = P(T_{\text{system}} > t)
\]

For Monte Carlo results, the tool can also compute:

- reliability confidence intervals
- MTTF and its confidence interval
- convergence behavior as simulation count increases
- contribution of critical paths to system failure

---

## Project Structure

```text
graph-based-reliability-analysis/
├── main.py
├── distributions.py
├── monte_carlo.py
├── critical_analysis.py
└── README.md
```

### File Descriptions

- `main.py` — GUI, workflow control, model management, and analysis execution
- `distributions.py` — supported lifetime distributions and reliability definitions
- `monte_carlo.py` — Monte Carlo simulation, convergence analysis, and component importance
- `critical_analysis.py` — plotting, validation summaries, hazard rate, and criticality analysis tools

---

## Installation

### Requirements

- Python 3.10 or higher
- NumPy
- SciPy
- SymPy
- Matplotlib
- PyQt6

Install dependencies:

```bash
pip install numpy scipy sympy matplotlib pyqt6
```

---

## Usage

Run the application:

```bash
python main.py
```

### Typical Workflow

1. Build a new graph-based model or load an existing JSON model
2. Select an analysis mode:
   - Static Analysis
   - Dynamic Reliability Analysis
   - Monte Carlo Simulation
3. Define component reliability parameters or lifetime distributions
4. Run the selected analysis
5. Inspect reliability curves, validation plots, and criticality outputs

---

## Example Outputs

The application can generate outputs such as:

- system reliability curve \(R(t)\)
- analytical vs Monte Carlo comparison plot
- Monte Carlo reliability curve with **95% CI**
- validation summary table
- lifetime histogram
- hazard rate plot
- critical interval summary
- path robustness comparison
- top-k critical paths
- path contribution plot
- component criticality / component importance charts

---

## Academic Context

This project was developed as an **undergraduate graduation project** in reliability engineering and system analysis.  
Its goal is to provide a flexible and interpretable framework for analyzing complex systems with heterogeneous component behaviors.

---

## Future Work

Possible future extensions include:

- rare-event acceleration methods
- Bayesian reliability modeling
- larger-scale system benchmarks
- optimization-based reliability improvement
- uncertainty quantification
- graph neural network-assisted reliability interpretation

---

## Author

**Selin Ayhan**  
Computer Engineering  
Ankara University

---

## License

This project is intended for **academic and educational use**.
