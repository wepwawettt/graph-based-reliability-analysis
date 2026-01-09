
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

---

## 🚀 Key Features

- Interactive graph-based system modeling  
- Automatic extraction of **minimal path sets**  
- Static and dynamic reliability analysis  
- Time-dependent reliability \( R(t) \)  
- Multiple lifetime distributions  
- Monte Carlo simulation  
- Common Cause Failure (CCF) modeling  
- Component Criticality Index (CCI)  
- Sensitivity analysis  
- PyQt6 GUI  

---

## 🧠 Methodology

1. Graph-based system modeling  
2. Minimal path extraction  
3. Analytical reliability (inclusion–exclusion)  
4. Time-dependent reliability computation  
5. Monte Carlo validation  
6. Criticality and robustness analysis  

---

## 📊 Monte Carlo Simulation

Empirical reliability function:

\[
R_{MC}(t) = P(T_{system} > t)
\]

---

## 🖥️ Project Structure

```text
├── main.py
├── distributions.py
├── monte_carlo.py
├── critical_analysis.py
├── README.md
```

---

## ⚙️ Installation

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

## ▶️ Usage

Run the application:

```bash
python main.py
```

Workflow:
1. Build or load a model  
2. Select analysis mode  
3. Run analysis  
4. Visualize results  

---

## 📈 Example Outputs

- System reliability curve \( R(t) \)  
- Monte Carlo reliability curve  
- Lifetime histogram  
- Critical interval plots  
- Component Criticality Index  

---

## 🔬 Academic Context

This project was developed as an **undergraduate graduation project** and follows standard reliability engineering methodologies.

---

## 🔮 Future Work

- Bayesian reliability  
- Graph Neural Networks  
- Large-scale optimization  
- Uncertainty quantification  

---

## 👩‍💻 Author

**Selin Ayhan**  
Computer Engineering  
Ankara University  

---

## 📄 License

For academic and educational use only.
````


