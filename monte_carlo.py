from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QPushButton, QLabel, QListWidget, QComboBox, QLineEdit,
    QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsTextItem,
    QGraphicsLineItem, QMessageBox, QFrame, QGraphicsItem, QGridLayout,
    QTextEdit, QMenu, QScrollArea , QDialog ,QDialogButtonBox ,QInputDialog, QTabWidget,QCheckBox,QGroupBox,QDoubleSpinBox,QSpinBox
    # QTabWidget buradan kaldırıldı
)
import numpy as np
from distributions import DISTRIBUTIONS
from PyQt6.QtWidgets import QFileDialog
from scipy.stats import lognorm, gamma

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure 
from PyQt6.QtGui import (
    QColor, QBrush, QPen, QFont, QPainter, QPainterPath, QPolygonF
)
from PyQt6.QtCore import Qt, QPointF, QObject
from PyQt6.QtWidgets import QInputDialog ,QInputDialog
from critical_analysis import plot_critical_intervals

# monte_carlo.py
import numpy as np

def run_monte_carlo(components, component_paths, N, t_max, ccf=None, n_t=100):
    """
    Monte Carlo simulation
    Returns:
      T_sys : array of system lifetimes
      t_vals: time grid
      R_mc  : Monte Carlo R(t)
    """
    T_sys = np.zeros(N, dtype=float)

    for i in range(N):
        lifetimes = {}
        T_ccf = None

        if ccf:
            beta, lambdas = ccf
            if beta > 0 and lambdas:
                lambda_avg = np.mean(lambdas)
                T_ccf = np.random.exponential(1 / (beta * lambda_avg))

        for cname, d in components.items():
            if d["dist"] == "static":
                lt_ind = 1e20
            else:
                conf = DISTRIBUTIONS[d["dist"]]
                lt_ind = conf["sample"](d["params"])

            lifetimes[cname] = min(lt_ind, T_ccf) if T_ccf else lt_ind

        path_fail_times = [
            min(lifetimes[c] for c in p)
            for p in component_paths if p
        ]

        T_sys[i] = max(path_fail_times) if path_fail_times else 1e20

    # 🔹 Monte Carlo R(t)
    t_vals = np.linspace(0, t_max, n_t)
    R_mc = np.array([np.mean(T_sys > t) for t in t_vals])

    return T_sys, t_vals, R_mc
