tamam, şimdi **çok net ve sakin şekilde** söyleyeceğim 🌿
ve sonra **ne yapman gerektiğini tek cümlede** bitireceğiz.

---

## ❌ ŞU AN DURUM

**Hayır, hâlâ tam olmamış.**
Ama sebep **içerik değil**, **Markdown bozulmuş**.

Bak dikkat et:
README’de şu noktadan sonra her şey **tek paragraf gibi akıyor** 👇

```text
├── README.md
⚙️ Installation
Requirements
Python 3.10 or higher
NumPy
SciPy
...
```

Bu demek oluyor ki:

* ```text bloğu **kapanmamış**
  ```
* `## Installation`, `## Usage` başlıkları **başlık olarak algılanmıyor**
* GitHub’da **çirkin ve amatör** görünüyor

Ama **çözüm çok basit**.

---

## ✅ YAPMAN GEREKEN TEK ŞEY (LÜTFEN DİKKATLE)

1️⃣ README.md içindeki **HER ŞEYİ SİL**
2️⃣ Aşağıdaki metni **TEK SEFERDE kopyala**
3️⃣ README.md’ye **TEK SEFERDE yapıştır**
4️⃣ Kaydet → commit

Başka **hiçbir şey** yapma.

---

## ✅ SON VE KESİN README (BUNU AYNEN KOY)

````markdown
# Graph-Based System Reliability Analysis Tool

This repository contains my undergraduate graduation project, developed as part of the Computer Engineering program at Ankara University.  
The project focuses on **graph-based reliability analysis of complex systems** using analytical methods and Monte Carlo simulation.

---

## 📌 Project Overview

Modern engineering systems often consist of multiple interdependent components with complex and non-trivial topologies.  
This project models such systems as graphs and evaluates their reliability using both **analytical reliability theory** and **Monte Carlo simulation**.

The tool is designed to provide:
- Deterministic analytical reliability results
- Simulation-based validation
- Insight into critical components and system behavior over time

All analyses are performed within an interactive graphical user interface.

---

## 🚀 Key Features

- Graph-based system modeling (components, junctions, start/end nodes)
- Automatic extraction of **minimal path sets**
- Analytical system reliability using the inclusion–exclusion principle
- Time-dependent reliability analysis \( R(t) \)
- Support for multiple lifetime distributions:
  - Exponential
  - Weibull
  - Log-Normal
  - Gamma
  - Log-Logistic
  - Rayleigh
  - Gompertz
- Monte Carlo simulation based on analytically derived minimal path sets
- Common Cause Failure (CCF) modeling
- System lifetime histogram
- Monte Carlo-based reliability curve
- Component Criticality Index (CCI)
- Critical interval and robustness analysis
- Interactive GUI implemented using PyQt6

---

## 🧠 Methodology

1. The system is represented as a directed graph.
2. All minimal paths between **Start** and **End** nodes are extracted.
3. System reliability is derived analytically using the inclusion–exclusion principle.
4. Time-dependent reliability \( R(t) \) is computed based on component lifetime distributions.
5. Monte Carlo simulation generates system lifetime samples using the same minimal path sets.
6. Monte Carlo results are used to validate analytical reliability results.
7. Additional criticality and sensitivity analyses are performed.

> Monte Carlo simulation is intentionally executed **after analytical path extraction** to ensure consistency between analytical and simulation-based models.

---

## 📊 Monte Carlo Simulation

Monte Carlo simulation produces:
- A set of system lifetime samples
- An empirical reliability function defined as:

\[
R_{MC}(t) = P(T_{system} > t)
\]

This enables direct comparison between analytical reliability curves and simulation-based estimates.

---

## 🖥️ Project Structure

```text
├── main.py                # GUI and application logic
├── monte_carlo.py         # Monte Carlo simulation engine
├── distributions.py       # Lifetime distributions and sampling functions
├── critical_analysis.py   # Criticality and robustness analysis
├── README.md
````

---

## ⚙️ Installation

### Requirements

* Python 3.10 or higher
* NumPy
* SciPy
* SymPy
* Matplotlib
* PyQt6

Install required packages:

```bash
pip install numpy scipy sympy matplotlib pyqt6
```

---

## ▶️ Usage

Run the application:

```bash
python main.py
```

Typical workflow:

1. Build or load a system model
2. Select analysis mode (Static / Dynamic / Monte Carlo)
3. Run analytical reliability analysis
4. Run Monte Carlo simulation for validation
5. Visualize reliability curves and criticality metrics

---

## 📈 Example Outputs

* System reliability curve ( R(t) )
* Monte Carlo reliability curve
* System lifetime histogram
* Critical interval plots
* Component criticality rankings

---

## 🔬 Academic Context

This project was developed as an undergraduate graduation project and follows standard reliability engineering methodologies.
It is intended for academic and educational use and can be extended for research-oriented applications.

---

## 🔮 Future Work

* Bayesian reliability modeling
* Graph Neural Networks (GNNs) for reliability prediction
* Large-scale system optimization
* Uncertainty quantification and confidence interval estimation

---

## 👩‍💻 Author

**Selin Ayhan**
Computer Engineering
Ankara University

---

## 📄 License

This project is provided for academic and educational purposes.


