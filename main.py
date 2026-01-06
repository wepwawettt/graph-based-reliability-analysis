"""
Graph-Based System Reliability Analysis Tool
Graduation Project – Ankara University
Author: Selin Ayhan
"""


import sys
import math
import numpy as np
from itertools import combinations
import os
import json
import time
import re
import copy
from critical_analysis import (
    plot_critical_intervals,
    plot_path_robustness,
    plot_critical_slopes,
    critical_summary_table,
    plot_component_criticality,
    plot_sensitivity_tornado,
    plot_true_path_robustness,
    plot_path_rt_curves,
    plot_survival_and_cdf,
    plot_system_rt_comparison
)

# --- GEREKLİ KÜTÜPHANELER ---
try:
    import sympy
    # LogNormal hatasını çözmek için erf, erfc, sqrt eklendi
    from sympy import erf, erfc, sqrt 
    from sympy import symbols, expand, exp, log, Symbol
except ImportError:
    print("HATA: 'sympy' kütüphanesi bulunamadı.")
    print("Lütfen 'py -3.11 -m pip install sympy' komutunu çalıştırdığınızdan emin olun.")
    sys.exit()
try:
    from scipy.stats import norm
    # LogNormal hatasını çözmek için 'scipy_erfc' eklendi
    from scipy.special import erfc as scipy_erfc 
except ImportError:
    print("UYARI: 'scipy' kütüphanesi bulunamadı. 'pip install scipy'")
    class norm:
        @staticmethod
        def cdf(*args, **kwargs): return 0.5
    # scipy yoksa, 'scipy_erfc'yi math kütüphanesinden almayı dene (alternatif)
    try:
        from math import erfc as scipy_erfc
    except ImportError:
        def scipy_erfc(*args, **kwargs): return 0.0 # Hiçbiri yoksa
# --- GEREKLİ KÜTÜPHANELER SONU ---

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QPushButton, QLabel, QListWidget, QComboBox, QLineEdit,
    QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsTextItem,
    QGraphicsLineItem, QMessageBox, QFrame, QGraphicsItem, QGridLayout,
    QTextEdit, QMenu, QScrollArea , QDialog ,QDialogButtonBox ,QInputDialog, QTabWidget,QCheckBox,QGroupBox,QDoubleSpinBox,QSpinBox
    # QTabWidget buradan kaldırıldı
)
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

from distributions import DISTRIBUTIONS

from monte_carlo import run_monte_carlo



# --- GÖRSEL SINIFLAR (Düğümler ve Düz Çizgiler) ---
class DraggableNode(QGraphicsEllipseItem):
    """ Sürüklenen 'Bileşen' (Node) """
    def __init__(self, *args, node_name=None, main_window=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.node_name = node_name; self.main_window = main_window
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        self.setAcceptHoverEvents(True)
        
    def mousePressEvent(self, event):
        self.setZValue(10); super().mousePressEvent(event)
    def mouseReleaseEvent(self, event):
        self.setZValue(0); super().mouseReleaseEvent(event)
        if self.main_window and self.node_name:
            self.main_window.update_node_position(self.node_name, self.pos())
    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)
        if self.main_window and self.node_name:
            self.main_window.update_node_position(self.node_name, self.pos())
    def hoverEnterEvent(self, event):
        QApplication.setOverrideCursor(Qt.CursorShape.PointingHandCursor); super().hoverEnterEvent(event)
    def hoverLeaveEvent(self, event):
        QApplication.restoreOverrideCursor(); super().hoverLeaveEvent(event)
    def contextMenuEvent(self, event):
        menu = QMenu()
        delete_action = menu.addAction("Bileşeni Sil")

        action = menu.exec(event.screenPos())
        if action == delete_action:
            self.main_window.remove_component(self.node_name)


class ClickableLineItem(QGraphicsLineItem):
    """ 
    Sağ tıklanabilen ve "Kavşak Ekle" menüsü çıkaran
    düz çizgi (tel).
    """
    def __init__(self, node1, node2, main_window, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.node1 = node1 # Başlangıç düğümü adı
        self.node2 = node2 # Bitiş düğümü adı
        self.main_window = main_window
        self.setPen(QPen(QColor("gray"), 2, Qt.PenStyle.SolidLine))
        self.setZValue(-1) # Düğümlerin arkasında kalması için
        self.setAcceptHoverEvents(True)

    def hoverEnterEvent(self, event):
        self.setPen(QPen(QColor("darkcyan"), 4, Qt.PenStyle.SolidLine))
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self.setPen(QPen(QColor("gray"), 2, Qt.PenStyle.SolidLine))
        super().hoverLeaveEvent(event)

    def contextMenuEvent(self, event):
        """ Sağ tıklandığında menü çıkar """
        menu = QMenu()
        split_action = menu.addAction(f"'{self.node1}' <-> '{self.node2}' arasına Kavşak Ekle")
        
        selected_action = menu.exec(event.screenPos())
        
        if selected_action == split_action:
            self.main_window.split_line(self.node1, self.node2, event.scenePos())

# --- GRAFİK PENCERESİ ---
class PlotWindow(QMainWindow):
    def __init__(self, t_values, plot_data, mttf=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Zamana Bağlı Güvenilirlik Grafiği R(t)")
        self.setGeometry(200, 200, 900, 600)

        self.figure = Figure(figsize=(8, 5), dpi=100, layout='constrained')
        self.canvas = FigureCanvas(self.figure)
        self.setCentralWidget(self.canvas)

        ax = self.figure.add_subplot(111)

        max_t = t_values[-1]

        # --- Eğrileri çiz ---
        for label, r_values in plot_data.items():
            label_with_value = f"{label} [R({max_t:.0f}) = {r_values[-1]:.4f}]"

            if label == "Sistem":
                ax.plot(
                    t_values, r_values,
                    label=label_with_value,
                    linewidth=2.8,        # daha kalın
                    color="navy",
                    zorder=10,
                )

            elif label.startswith("Yol_"):
                ax.plot(
                    t_values, r_values,
                    label=label_with_value,
                    linewidth=1.0,        # ince
                    linestyle="--",
                    alpha=0.45,           # soluk renk
                    zorder=1,
                )

            else:
                ax.plot(
                    t_values, r_values,
                    label=label_with_value,
                    linewidth=1.6,        # bileşenler orta kalınlık
                    linestyle="-.",
                    alpha=0.9,
                    zorder=5,
                )


        ax.set_xlabel("Zaman (t)")
        ax.set_ylabel("Güvenilirlik R(t)")
        ax.set_title("Sistem ve Bileşenlerin Güvenilirlik Analizi R(t)")
        ax.set_ylim(0, 1.05)
        ax.set_xlim(0, max_t)
        ax.grid(True, linestyle='--', linewidth=0.5)

        # ============================================
        #  SISTEM EĞRİSİNİ AL
        # ============================================
        if "Sistem" in plot_data:
            system_r = plot_data["Sistem"]
        else:
            system_r = list(plot_data.values())[0]

        # ============================================
        #  t_90 VE t_10 NOKTALARINI BUL
        # ============================================
        t_90, t_10 = None, None

        for i in range(len(t_values)):
            if t_90 is None and system_r[i] <= 0.9:
                t_90 = t_values[i]
            if t_10 is None and system_r[i] <= 0.1:
                t_10 = t_values[i]
                break

        # ============================================
        #  KRITIK BÖLGE (0.1 – 0.9)
        # ============================================
        ax.axhline(0.9, color='green', linestyle='--', alpha=0.6)
        ax.axhline(0.1, color='red', linestyle='--', alpha=0.6)

        ax.fill_between(t_values, 0.1, 0.9,
                        color='yellow', alpha=0.15,
                        label="Kritik Bölge (0.1 - 0.9)")

        # ============================================
        #  t_90 DİKEY ÇİZGİ
        # ============================================
        if t_90:
            ax.axvline(t_90, color='green', linestyle=':', linewidth=2)
            ax.text(t_90, 0.93,
                    f"t₉₀ = {t_90:.1f}",
                    rotation=90,
                    color='green', fontsize=10,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="green"))

        # ============================================
        #  t_10 DİKEY ÇİZGİ
        # ============================================
        if t_10:
            ax.axvline(t_10, color='red', linestyle=':', linewidth=2)
            ax.text(t_10, 0.15,
                    f"t₁₀ = {t_10:.1f}",
                    rotation=90,
                    color='red', fontsize=10,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="red"))

        # ============================================
        #  MTTF DİKEY ÇİZGİ
        # ============================================
        if mttf is not None:
            ax.axvline(mttf, color='purple', linestyle='-.', linewidth=2)
            ax.text(mttf, 0.50,
                    f"MTTF = {mttf:.1f}",
                    rotation=90,
                    color='purple', fontsize=10,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="purple"))

        # ============================================
        #  BİLGİ KUTUSU (t90, t10, MTTF)
        # ============================================
        info = []
        if t_90: info.append(f"t₉₀ (R=0.9): {t_90:.1f}")
        if t_10: info.append(f"t₁₀ (R=0.1): {t_10:.1f}")
        if mttf is not None: info.append(f"MTTF: {mttf:.1f}")

        ax.text(0.98, 0.02,
                "\n".join(info),
                transform=ax.transAxes,
                fontsize=10,
                ha='right', va='bottom',
                bbox=dict(facecolor="white", alpha=0.8))

        ax.legend()
        self.canvas.draw()



class HistogramWindow(QMainWindow):
    def __init__(self, data, bins=10, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sistem Ömrü Histogramı (Analitik + Simülasyon)")
        self.setGeometry(200, 200, 800, 500)

        self.figure = Figure(figsize=(7, 4), dpi=100, layout='constrained')
        self.canvas = FigureCanvas(self.figure)
        self.setCentralWidget(self.canvas)

        ax = self.figure.add_subplot(111)

        # === MATLAB ile birebir aynı histogram ===
        # 1) Bin genişliğini 5 yapıyoruz (0–5–10–15...)
        # === Kitap stilinde histogram ===
        # === Kitap Formatında Histogram (Şekil 5.15 gibi) ===

# 1) Sistemin maksimum arıza zamanı t_max (ör. 90000)
        max_val = max(data)

        # 2) Kitaptaki aralığa ölçekle: 0–50 arası olsun
        #    Bu şekilde kitap histogramına birebir benzer.
        scale_factor = max_val / 50      # istediğin eksen uzunluğu kadar bölebilirsin
        data_scaled = data / scale_factor

        # 3) Kitaptaki gibi bin genişliği = 5
        bin_edges = np.arange(0, 55, 5)   # 0,5,10,...,50

        ax.hist(
            data_scaled,
            bins=bin_edges,
            color="gray",
            edgecolor="black",
            linewidth=0.8
        )

        # 4) X ekseni kitap formatında olsun
        ax.set_xticks(np.arange(0, 55, 5))

        ax.set_xlabel("Sistem Arıza Zamanı (ölçeklenmiş) [Kitap Formatı]")
        ax.set_ylabel("Frekans")
        ax.set_title("Simüle Edilmiş Sistem Ömrü Histogramı")

# --- FORMÜL PENCERESİ (v13 - TEK SAYFA KAYDIRILABİLİR GÖRÜNÜM) ---
class FormulaWindow(QMainWindow):
    def __init__(self, formula_latex_list, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Üretilen Analitik Formüller (Açıklamalı)")
        self.setGeometry(250, 250, 1000, 700) 

        self.scroll_area = QScrollArea()
        self.setCentralWidget(self.scroll_area)
        
        # Dikey kaydırma için False olmalı (v7 düzeltmesi)
        self.scroll_area.setWidgetResizable(False) 

        # Toplam satır sayısını ve en uzun satırı bul
        all_lines = []
        for block in formula_latex_list:
            all_lines.extend(block.split('\n'))
        
        all_lines_cleaned = [line.replace('$', '') for line in all_lines]
        total_lines = len(all_lines)
        max_len = 0
        if all_lines_cleaned:
             max_len = max(len(line) for line in all_lines_cleaned)

        fig_width_inches = max(12, max_len * 0.1) 
        fig_height_inches = max(8, total_lines * 0.35) # Satır aralığı biraz artırıldı
        
        self.formula_figure = Figure(figsize=(fig_width_inches, fig_height_inches), dpi=100)
        self.formula_figure.text(0.01, 0.99, "", usetex=False)
        self.formula_canvas = FigureCanvas(self.formula_figure)
        
        self.scroll_area.setWidget(self.formula_canvas)

        self.formula_ax = self.formula_figure.add_subplot(111)
        self.update_formula_display(formula_latex_list)
    def wrap_math(self, s):
        return r"$" + s + r"$"

    def update_formula_display(self, formula_latex_list):
        self.formula_ax.clear()
        self.formula_ax.axis('off')
        formula_latex_list = [f.replace("$", "") for f in formula_latex_list]
        y = 0.98
        dy = 0.06   # satır aralığı

        # --- 1. Başlık ---
        self.formula_ax.text(0.01, y, "--- 1. Geniş Formül (Yol Bazlı) ---",
                            fontsize=12, ha="left", va="top")
        y -= dy

        # --- Yol bazlı formül ---
        for line in formula_latex_list[0].split("\n"):
            line = line.strip()
            if not line:
                continue   # ⬅️ BOŞ SATIRI ATLAR
            self.formula_ax.text(0.03, y, rf"${line}$", fontsize=12)
            y -= dy


        # --- 2. Başlık ---
        self.formula_ax.text(0.01, y, "--- 2. Sade Formül (Bileşen Bazlı) ---",
                            fontsize=12, ha="left", va="top")
        y -= dy

        # --- Bileşen formülü ---
        for line in formula_latex_list[1].split("\n"):
            line = line.strip()
            if not line:
                continue
            self.formula_ax.text(0.03, y, rf"${line}$",
                                fontsize=12, ha="left", va="top")
            y -= dy

            

        # --- 3. Başlık ---
        self.formula_ax.text(0.01, y, "--- 3. Zamana Bağlı Formül (R(t)) ---",
                            fontsize=12, ha="left", va="top")
        y -= dy

        # --- R(t) formülü ---
        for line in formula_latex_list[2].split("\n"):
            line = line.strip()
            if not line:
                continue
            self.formula_ax.text(0.03, y, rf"${line}$",
                                fontsize=12, ha="left", va="top")
            y -= dy

        self.formula_canvas.draw()



# --- TEK SAYFA FORMÜL PENCERESİ SONU ---


# --- ANA PENCERE (Hibrit Model: Bileşen + Otomatik Kavşak) ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Analitik Güvenilirlik Çözücü (v13 - Tek Sayfa Formül)")
        self.setGeometry(100, 100, 1300, 850)
        

        # --- Varsayılan analiz modu ---
        self.analysis_mode = "static"

        # --- Veri Modeli ---
        self.components = {}      # Bozulabilenler (a1, a5...)
        self.junctions = set()    # Bozulmaz kavşaklar (j1, j2...)
        self.junction_count = 1   # j1, j2... diye isim vermek için
        self.graph = {}           # Çift yönlü mantık için
        self.connections = []

        self.node_positions = {}
        self.node_items = {}
        self.edge_items = {}

        self.plot_window = None
        self.formula_window = None
        self.formula_latex = None  # Formüller LaTeX metin blokları

        # === HESAPLAMA WIDGET'LARI (sağ panelde kullanılacak) ===
        
        self.run_button = QPushButton("FORMÜL ÜRET & HESAPLA")
        self.run_button.setMinimumHeight(36)
        self.run_button.setStyleSheet("background-color: lightgreen; font-weight: bold;")
        self.run_button.clicked.connect(self.run_main_button)

        self.result_label = QLabel("Sistem Güvenirliği: -")
        self.result_label.setFont(QFont("Arial", 12, QFont.Weight.Bold))

        self.show_formula_button = QPushButton("Üretilen Formülleri Göster")
        self.show_formula_button.clicked.connect(self.show_formula_window)
        self.show_formula_button.setEnabled(True)


        # === MERKEZ WIDGET & ANA LAYOUT ===
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # === SOL PANEL: MODEL KURMA PANELİ ===
        # === SOL PANEL (ScrollArea'lı) ===
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(350)
        scroll.setMaximumWidth(350)

        model_builder_panel = QWidget()
        model_layout = QVBoxLayout(model_builder_panel)
        scroll.setWidget(model_builder_panel)

        main_layout.addWidget(scroll)


        # 1. Adım: Analiz Modu
        model_layout.addWidget(QLabel("1. Adım: Analiz Modunu Seçin"))
        self.mode_selector = QComboBox()
        self.mode_selector.addItems([
            "Statik Analiz (R)",
            "Dinamik Analiz R(t)",
            "Monte Carlo Simülasyonu"
        ])
        self.mode_selector.setMinimumHeight(32)
        self.mode_selector.currentTextChanged.connect(self.on_mode_changed)
        model_layout.addWidget(self.mode_selector)

        # 2. Adım: Bileşen Ekleme
        self.comp_name_input = QLineEdit()
        self.comp_name_input.setPlaceholderText("Bileşen Adı (örn: a1, a2...)")
        self.comp_name_input.setMinimumHeight(32)
        model_layout.addWidget(self.comp_name_input)

        # --- STATİK GİRİŞLER ---
        self.static_inputs_widget = QWidget()
        static_layout = QVBoxLayout(self.static_inputs_widget)
        static_layout.setContentsMargins(0, 0, 0, 0)

        self.comp_reli_input = QLineEdit()
        self.comp_reli_input.setPlaceholderText("Güvenirlik (R) (örn: 0.99)")
        self.comp_reli_input.setMinimumHeight(32)
        static_layout.addWidget(self.comp_reli_input)

        model_layout.addWidget(self.static_inputs_widget)

        # --- DİNAMİK GİRİŞLER ---
        self.dynamic_inputs_widget = QWidget()
        dynamic_layout = QGridLayout(self.dynamic_inputs_widget)
        dynamic_layout.setContentsMargins(0, 5, 0, 0)
        self.param_inputs = []
        self.param_layout = QVBoxLayout()
        dynamic_layout.addLayout(self.param_layout, 3, 0, 1, 2)


        self.dist_label = QLabel("Dağılım Tipi:")
        self.dist_selector = QComboBox()
        self.dist_selector.addItems(["Exponential", "Weibull", "Log-Normal",
                             "Gamma", "Log-Logistic", "Gompertz", "Rayleigh"])

        self.dist_selector.setMinimumHeight(32)
        self.dist_selector.currentTextChanged.connect(self.on_dist_changed)

        dynamic_layout.addWidget(self.dist_label, 0, 0)
        dynamic_layout.addWidget(self.dist_selector, 0, 1)

        self.param1_label = QLabel("Ölçek (Scale/θ):")
        self.param1_input = QLineEdit()
        self.param1_input.setPlaceholderText("örn: 10000 (1/λ)")
        self.param1_input.setMinimumHeight(32)
        dynamic_layout.addWidget(self.param1_label, 1, 0)
        dynamic_layout.addWidget(self.param1_input, 1, 1)

        self.param2_label = QLabel("Şekil (k/β):")
        self.param2_input = QLineEdit()
        self.param2_input.setPlaceholderText("örn: 1.5")
        self.param2_input.setMinimumHeight(32)
        dynamic_layout.addWidget(self.param2_label, 2, 0)
        dynamic_layout.addWidget(self.param2_input, 2, 1)
        self.param1_input.hide()
        self.param2_input.hide()
        self.param1_label.hide()
        self.param2_label.hide()

        model_layout.addWidget(self.dynamic_inputs_widget)

        # Dinamik başlangıç durumu
        self.dynamic_inputs_widget.hide()
        self.on_dist_changed()
        

        # Bileşen ekleme butonları / liste
        self.add_comp_button = QPushButton("Bileşen Ekle")
        self.add_comp_button.clicked.connect(self.add_component)
        model_layout.addWidget(self.add_comp_button)

        self.comp_list_widget = QListWidget()
        model_layout.addWidget(self.comp_list_widget)

        
        self.update_comp_button = QPushButton("Seçili Bileşeni Güncelle")
        model_layout.addWidget(self.update_comp_button)
        self.update_comp_button.setEnabled(False)
        self.comp_list_widget.currentRowChanged.connect(
        self.on_component_selected
    )
        self.update_comp_button.clicked.connect(self.update_selected_component)

        # 3. Adım: Bağlantı Ekleme (Teller)
        model_layout.addWidget(self.create_separator("3. Adım: Bağlantıları (Telleri) Ekle"))
        model_layout.addWidget(QLabel("NOT: Bir 'tele' bağlanmak için, sahnedeki tele sağ tıklayın."))

        self.from_node_selector = QComboBox()
        self.to_node_selector = QComboBox()
        self.from_node_selector.setMinimumHeight(32)
        self.to_node_selector.setMinimumHeight(32)

        self.add_connection_button = QPushButton("Bağlantı Ekle (Düğüm 1 <-> Düğüm 2)")
        self.add_connection_button.clicked.connect(self.add_connection)

        model_layout.addWidget(QLabel("Düğüm 1:"))
        model_layout.addWidget(self.from_node_selector)
        model_layout.addWidget(QLabel("Düğüm 2:"))
        model_layout.addWidget(self.to_node_selector)
        model_layout.addWidget(self.add_connection_button)

        self.conn_list_widget = QListWidget()
        model_layout.addWidget(self.conn_list_widget)

        # 5. Adım: Formül Göstergesi (başlık sadece)
        model_layout.addWidget(self.create_separator("5. Adım: Üretilen Formüller"))

        # Model kaydet / yükle butonları
        self.save_model_button = QPushButton("Modeli Kaydet")
        
        model_layout.addWidget(self.save_model_button)

        self.load_model_button = QPushButton("Model Yükle")
        self.load_model_button.clicked.connect(self.load_model_dialog)
        self.save_model_button.clicked.connect(self.save_model_dialog)

        model_layout.addWidget(self.load_model_button)

        model_layout.addStretch()

        # === ORTA PANEL: SEKME BAZLI GRAFİK SAHNELER ===
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget, 10)
        self.tab_widget.currentChanged.connect(self.on_tab_changed)
        # === BAŞLANGIÇTA BOŞ (MANUEL) MODEL SEKME ===
        self.create_empty_model_tab("Yeni Model")


        # === SAĞ PANEL: HESAPLAMA & SONUÇLAR ===
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(10, 10, 10, 10)

        right_layout.addWidget(self.create_separator("4. Adım: Hesaplama ve Sonuçlar"))
       # === MONTE CARLO DENEME SAYISI ===
        mc_box = QGroupBox("Monte Carlo Ayarları")
        mc_layout = QVBoxLayout()

        self.mc_spinbox = QSpinBox()
        self.mc_spinbox.setRange(100, 1_000_000)
        self.mc_spinbox.setSingleStep(1000)
        self.mc_spinbox.setValue(20000)

        mc_layout.addWidget(QLabel("Simülasyon Sayısı (N):"))
        mc_layout.addWidget(self.mc_spinbox)

        mc_box.setLayout(mc_layout)
        right_layout.addWidget(mc_box)

      
        # === CCF (Common Cause Failure) ===
        self.ccf_checkbox = QPushButton("CCF Kullan (β)")
        self.ccf_checkbox.setCheckable(True)
        self.ccf_checkbox.setChecked(False)

        self.ccf_beta_input = QLineEdit("0.0")
        self.ccf_beta_input.setPlaceholderText("β (0–1)")
        self.ccf_beta_input.setEnabled(False)

        self.ccf_checkbox.toggled.connect(
            lambda checked: self.ccf_beta_input.setEnabled(checked)
        )

        right_layout.addWidget(self.ccf_checkbox)
        right_layout.addWidget(self.ccf_beta_input)
        # === ZAMAN ARALIĞI (T) ===
        time_box = QGroupBox("Zaman Ayarları")
        time_layout = QVBoxLayout()

        self.t_max_input = QDoubleSpinBox()
        self.t_max_input.setRange(1, 1e7)
        self.t_max_input.setValue(1000)
        self.t_max_input.setSuffix(" s")
        self.t_max_input.setDecimals(0)

        time_layout.addWidget(QLabel("Maksimum Zaman (t_max):"))
        time_layout.addWidget(self.t_max_input)

        time_box.setLayout(time_layout)
        right_layout.addWidget(time_box)

        right_layout.addWidget(self.run_button)
        right_layout.addWidget(self.result_label)
        right_layout.addWidget(self.show_formula_button)
        right_layout.addStretch()
        self.critical_button = QPushButton("Critical Analysis (Multi-Model)")
        self.critical_button.setMinimumHeight(36)
        self.critical_button.clicked.connect(self.run_critical_analysis)
        right_layout.addWidget(self.critical_button)

        main_layout.addWidget(right_panel)

        # Başlangıçta mod ayarını bir kez uygula
        self.on_mode_changed("Statik Analiz (R)")
        self.setup_connection_context_menu()
        self.comp_list_widget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.comp_list_widget.customContextMenuRequested.connect(
            self.show_component_context_menu
        )
    def show_component_context_menu(self, pos):
        item = self.comp_list_widget.itemAt(pos)
        if not item:
            return

        comp_name = item.text()

        menu = QMenu()
        delete_action = menu.addAction("Bileşeni Sil")

        action = menu.exec(self.comp_list_widget.mapToGlobal(pos))
        if action == delete_action:
            self.remove_component(comp_name)

    def on_tab_changed(self, index):
        tab = self.tab_widget.widget(index)
        if not tab or not hasattr(tab, "model_state"):
            return

        state = tab.model_state

        # 🔁 MODEL STATE BAĞLA
        self.components      = state["components"]
        self.junctions       = state["junctions"]
        self.connections     = state["connections"]
        self.graph           = state["graph"]
        self.node_positions  = state["node_positions"]
        self.node_items      = state["node_items"]
        self.edge_items      = state["edge_items"]

        # 🎨 SCENE / VIEW
        view = tab.findChild(QGraphicsView)
        if view:
            self.view = view
            self.scene = view.scene()

        # 🔄 SOL PANEL
        self.refresh_left_panel()


    def run_critical_analysis(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Critical Analysis – Model Seçimi")
        dialog.resize(300, 400)

        layout = QVBoxLayout(dialog)

        label = QLabel("Karşılaştırmak istediğiniz modelleri seçin:")
        layout.addWidget(label)

        checkboxes = []

        # === TABLARDAN MODEL LİSTESİ ===
        for i in range(self.tab_widget.count()):
            name = self.tab_widget.tabText(i)
            tab = self.tab_widget.widget(i)

            # sadece analiz yapılmış modeller
            if not hasattr(tab, "model_state"):
                continue
            if "analysis_results" not in tab.model_state:
                continue

            cb = QCheckBox(name)
            layout.addWidget(cb)
            checkboxes.append((cb, tab))

        # === BUTONLAR ===
        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("Analiz Et")
        cancel_btn = QPushButton("İptal")

        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

        cancel_btn.clicked.connect(dialog.reject)

        def on_ok():
            selected = [
                (cb.text(), tab)
                for cb, tab in checkboxes
                if cb.isChecked()
            ]

            if len(selected) < 2:
                QMessageBox.warning(
                    dialog,
                    "Yetersiz Seçim",
                    "Critical Analysis için en az 2 model seçmelisiniz."
                )
                return

            dialog.accept()
            self._run_critical_analysis_on_selected(selected)

        ok_btn.clicked.connect(on_ok)

        dialog.exec()

    def _run_critical_analysis_on_selected(self, selected_tabs):
        results_dict = {}
        paths_dict = {}
        path_rt_dict = {}

        for name, tab in selected_tabs:
            analysis = tab.model_state["analysis_results"]

            results_dict[name] = {
                "t": analysis["t"],
                "R": analysis["R"]
            }

            paths = tab.model_state.get("component_paths", [])
            path_rts = tab.model_state.get("path_rt", [])

            # sadece path'i olan modelleri ekle
            if paths and path_rts:
                paths_dict[name] = paths
                path_rt_dict[name] = path_rts

        # === MODEL SEVİYESİ ANALİZLER ===
        # ✅ Multi-model R(t) tek grafikte
        plot_system_rt_comparison(results_dict)

        intervals = plot_critical_intervals(results_dict)
        plot_critical_slopes(intervals)

        df = critical_summary_table(intervals)
        print("\n=== CRITICAL ANALYSIS SUMMARY ===")
        print(df)

        # === PATH SEVİYESİ ANALİZLER ===
        if paths_dict:
            plot_path_rt_curves(results_dict, paths_dict, path_rt_dict)
            plot_true_path_robustness(results_dict, paths_dict, path_rt_dict)
                    # === COMPONENT CRITICALITY (CCI) ===
            # === COMPONENT CRITICALITY (MODEL BAZLI) ===
            for model_name in results_dict:
                single_results = {
                    model_name: results_dict[model_name]
                }
                single_paths = {
                    model_name: paths_dict[model_name]
                }

                plot_component_criticality(single_results, single_paths)






    def create_empty_model_tab(self, tab_name):
        tab = QWidget()
        tab_layout = QVBoxLayout(tab)
        

        scene = QGraphicsScene()
        scene.setBackgroundBrush(QBrush(QColor(230, 230, 230)))

        view = QGraphicsView(scene)
        view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        view.setRenderHint(QPainter.RenderHint.Antialiasing)

        tab_layout.addWidget(view)
        tab.setLayout(tab_layout)

        self.tab_widget.addTab(tab, tab_name)
        tab.model_state = {
    "components": {},
    "junctions": set(),
    "connections": [],
    "graph": {},
    "node_positions": {},
    "node_items": {},
    "edge_items": {}
}
        # 👉 BU SEKMEYİ AKTİF YAP
        self.tab_widget.setCurrentWidget(tab)

        # 👉 SEKME STATE'İNİ AKTİF ET (ÇOK ÖNEMLİ)
        self.on_tab_changed(self.tab_widget.currentIndex())

        # 👉 ARTIK GÜVENLİ: Start / End çiz
        self.initialize_scene()


    # --- SİMÜLASYON MOTORU (Tüm düzeltmeler dahil) ---
    def run_analysis(self):
        mttf = None
        t_max = self.t_max_input.value()

        if not self.components or not self.graph:
            QMessageBox.warning(self, "Hata", "Önce bileşen ve bağlantı tanımlanmalı.")
            return

        self.run_button.setText("Hesaplanıyor...")
        QApplication.processEvents()

        try:
            # === 1. PATH SETS ===
            print("1. Tüm minimal yollar (path sets) bulunuyor...")
            paths_nodes = self._find_all_paths("Start", "End")

            if not paths_nodes:
                QMessageBox.critical(self, "Hata", "Start ile End arasında yol yok.")
                self.run_button.setText("FORMÜL ÜRET & HESAPLA")
                return

            component_paths = [
                frozenset([n for n in p if n in self.components])
                for p in paths_nodes
            ]
            component_paths = list(dict.fromkeys(component_paths))
            tab = self.tab_widget.currentWidget()
            tab.model_state["component_paths"] = component_paths

            print("  Bulunan yollar:", component_paths)

            # === 2. SEMBOLİK FORMÜLLER ===
            print("2. Sembolik formüller üretiliyor...")

            self.formula_latex = ["", "", ""]
            genis_formula_str_list = [
    r"\text{Sistemde " + str(len(component_paths)) + r" başarılı yol vardır:}"
]

            path_symbols = [symbols(f"P_{i+1}") for i in range(len(component_paths))]

            for i, pset in enumerate(component_paths):
                term = r" \cdot ".join([f"R_{{{c}}}" for c in sorted(list(pset))])
                genis_formula_str_list.append(f"  $P_{i+1}$ = {term}")

            p_union = r" \cup ".join([f"P_{i+1}" for i in range(len(component_paths))])
            genis_formula_str_list.append(f"$R_{{Sistem}} = P({p_union})$")

            wide_formula = 0
            for k in range(1, len(path_symbols) + 1):
                for comb in combinations(path_symbols, k):
                    term = 1
                    for p in comb:
                        term *= p
                    wide_formula += ((-1)**(k+1)) * term

            genis_formula_str_list.append(sympy.latex(wide_formula))
            self.formula_latex[0] = self.clean_latex("\n".join(genis_formula_str_list))

            # === 3. BİLEŞEN FORMÜLÜ ===
            comp_symbols = {c: symbols(f"R_{c}") for c in self.components}
            final_formula = 0

            for k in range(1, len(component_paths) + 1):
                for comb in combinations(range(len(component_paths)), k):
                    union_c = frozenset.union(*[component_paths[i] for i in comb])
                    term = 1
                    for cname in union_c:
                        term *= comp_symbols[cname]
                    final_formula += ((-1)**(k+1)) * term

            final_formula = expand(final_formula)
            self.formula_latex[1] = self.clean_latex(sympy.latex(final_formula))

            
            # === 4. DİNAMİK R(t) FORMÜLÜ ===
            t = symbols("t", positive=True)
            subs_dict_calc = {}
            comp_plot_data = {}

            for comp_name, data in self.components.items():
                if data["dist"] == "static":
                    rt = data["R"]
                else:
                    conf = DISTRIBUTIONS[data["dist"]]
                    if conf.get("R_sym") is not None:
                        rt = conf["R_sym"](t, data["params"])
                    elif "R_num" in conf:
                        # NUMERİK dağılımlar sembolik formüle girmez
                        rt = symbols(f"R_{comp_name}")

                    else:
                        rt = 1




                subs_dict_calc[symbols(f"R_{comp_name}")] = rt
                comp_plot_data[comp_name] = data

            # === 5. SAYISAL R(t) HESABI ===
            # === ZAMAN VEKÖRÜ (LOG-GÜVENLİ) ===
            
            t_safe = np.linspace(1e-6, t_max, 400)   # 0 YOK
            t_values = t_safe.copy()


            modules = ["numpy", {"exp": np.exp, "log": np.log, "sqrt": np.sqrt, "erfc": scipy_erfc}]
            final_rt_formula = final_formula.subs(subs_dict_calc)
            
            rt_latex = sympy.latex(final_rt_formula)
            self.formula_latex[2] = self.clean_latex(
                r"R_{\text{Sistem}}(t) = " + rt_latex
            )

            # 1) Önce sistem R(t)'yi 1 kabul et
            system_r = np.ones_like(t_safe, dtype=float)
            R_ccf_numeric = None

            if self.ccf_checkbox.isChecked():
                try:
                    beta = float(self.ccf_beta_input.text())
                    if not (0 <= beta <= 1):
                        raise ValueError
                except:
                    QMessageBox.warning(self, "Hata", "CCF β değeri 0–1 arasında olmalı.")
                    return

                lambdas = [
                    d["params"]["lambda"]
                    for d in self.components.values()
                    if d["dist"] == "Exponential"
                ]

                if lambdas:
                    lambda_avg = np.mean(lambdas)
                    R_ccf_numeric = np.exp(-lambda_avg * t_safe)


            # 2) Her minimal yol için numerik R(t) hesapla
            path_rts = []

            for pset in component_paths:
                rt_path = np.ones_like(t_safe, dtype=float)

                for cname in pset:
                    data = self.components[cname]

                    # --- 1) Bağımsız bileşen R(t) ---
                    if data["dist"] == "static":
                        rt_ind = np.ones_like(t_safe) * data["R"]
                    else:
                        conf = DISTRIBUTIONS[data["dist"]]
                        if conf.get("R_sym") is not None:
                            f = sympy.lambdify(t, conf["R_sym"](t, data["params"]), modules)
                            rt_ind = f(t_safe)
                        else:
                            rt_ind = conf["R_num"](t_safe, data["params"])

                    rt_ind = np.asarray(rt_ind, dtype=float)
                    rt_ind = np.nan_to_num(rt_ind, 0.0, 0.0, 0.0)

                    # --- 2) CCF uygula (ASLİ YERİ BURASI) ---
                    if R_ccf_numeric is not None:
                        rt_c = (1 - beta) * rt_ind + beta * R_ccf_numeric
                    else:
                        rt_c = rt_ind

                    rt_path *= rt_c

                path_rts.append(rt_path)

            # 3) Inclusion–Exclusion (numeric)
            system_r = np.zeros_like(t_safe)

            for k in range(1, len(path_rts) + 1):
                for comb in combinations(range(len(path_rts)), k):
                    prod = np.ones_like(t_safe)
                    for i in comb:
                        prod *= path_rts[i]
                    system_r += ((-1) ** (k + 1)) * prod

            system_r = np.clip(system_r, 0.0, 1.0)

            system_r = np.asarray(system_r, dtype=float)
            system_r = np.nan_to_num(system_r, nan=0.0, posinf=0.0, neginf=0.0)


            plot_data_final = {"Sistem": system_r}

            for cname, data in comp_plot_data.items():
                if data["dist"] == "static":
                    plot_data_final[cname] = np.ones_like(t_safe) * data["R"]
                else:
                    conf = DISTRIBUTIONS[data["dist"]]

                    if conf.get("R_sym") is not None:
                        f = sympy.lambdify(t, conf["R_sym"](t, data["params"]), modules)
                        plot_data_final[cname] = f(t_safe)
                        plot_data_final[cname] = np.asarray(plot_data_final[cname], dtype=float)
                        plot_data_final[cname] = np.nan_to_num(plot_data_final[cname], 0.0, 0.0, 0.0)

                    elif "R_num" in conf:
                        plot_data_final[cname] = conf["R_num"](t_safe, data["params"])
                        plot_data_final[cname] = np.asarray(plot_data_final[cname], dtype=float)
                        plot_data_final[cname] = np.nan_to_num(plot_data_final[cname], 0.0, 0.0, 0.0)

                    else:
                        plot_data_final[cname] = np.ones_like(t_safe)

                # === YOLLARIN R(t) EĞRİLERİ ===
            for i, pset in enumerate(component_paths):
                components_in_path = " → ".join(sorted(list(pset)))
                path_name = f"Yol_{i+1} ({components_in_path})"

                rt = np.ones_like(t_safe)
                for cname in pset:
                    data = self.components[cname]

                    if data["dist"] == "static":
                        rt_ind = np.ones_like(t_safe) * data["R"]
                    else:
                        conf = DISTRIBUTIONS[data["dist"]]
                        if conf.get("R_sym") is not None:
                            f = sympy.lambdify(t, conf["R_sym"](t, data["params"]), modules)
                            rt_ind = f(t_safe)
                        else:
                            rt_ind = conf["R_num"](t_safe, data["params"])

                    rt_ind = np.nan_to_num(rt_ind, 0.0, 0.0, 0.0)

                    if R_ccf_numeric is not None:
                        rt_c = (1 - beta) * rt_ind + beta * R_ccf_numeric
                    else:
                        rt_c = rt_ind

                    rt *= rt_c

                plot_data_final[path_name] = rt


            # 🔴 PATH R(t) LİSTESİNİ SAKLA (Critical Analysis için)
            tab = self.tab_widget.currentWidget()
            tab.model_state["component_paths"] = component_paths
            tab.model_state["path_rt"] = [
                plot_data_final[f"Yol_{i+1} ({' → '.join(sorted(list(pset)))})"]
                for i, pset in enumerate(component_paths)
            ]

            # === MTTF HESABI ===
            try:
                mttf = np.trapz(system_r, t_safe)
                if not np.isfinite(mttf):
                    mttf = 0.0

            except:
                mttf = None

            # === GRAFİK ===
            self.plot_window = PlotWindow(t_safe, plot_data_final, mttf=mttf)

            self.plot_window.show()
            # === ANALİZ SONUÇLARINI SAKLA (KRİTİK ANALİZ İÇİN) ===
            self.last_results = {
                "CurrentSystem": {
                    "t": t_safe,
                    "R": system_r
                }
            }


        except Exception as e:
            QMessageBox.critical(self, "Analiz Hatası", f"Hata: {e}")
            print("HATA:", e)

        self.run_button.setText("FORMÜL ÜRET & HESAPLA")

        mttf_text = f"{mttf:.2f}" if mttf is not None else "N/A"

        self.result_label.setText(
            f"Sistem Güvenirliği: R(t={t_max:.0f}) = {float(system_r[-1]):.6f},  MTTF ≈ {mttf_text}"
        )

                # === SENSITIVITY / TORNADO ANALYSIS ===
        try:
            mttf_base, sensitivity = self.run_sensitivity_analysis()
            if sensitivity:
                plot_sensitivity_tornado(mttf_base, sensitivity)
        except Exception as e:
            print("Sensitivity analysis hatası:", e)

        # === ANALİZ SONUCUNU AKTİF SEKMEYE KAYDET (ÇOK ÖNEMLİ) ===
        tab = self.tab_widget.currentWidget()
        tab.model_state["analysis_results"] = {
            "t": t_safe,
            "R": system_r
        }

        tab = self.tab_widget.currentWidget()
        tab.model_state["analysis_results"] = {
            "t": t_safe,
            "R": system_r
        }

    def on_dist_changed(self):
        # Eski parametre widgetlarını temizle
        for label, edit in getattr(self, "param_inputs", []):
            label.deleteLater()
            edit.deleteLater()
        self.param_inputs = []

        dist = self.dist_selector.currentText()
        config = DISTRIBUTIONS[dist]

        for prm in config["params"]:
            label = QLabel(prm["label"])
            edit = QLineEdit()
            edit.setMinimumHeight(28)

            self.param_layout.addWidget(label)
            self.param_layout.addWidget(edit)

            self.param_inputs.append((label, edit))



    def setup_connection_context_menu(self):
        self.conn_list_widget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.conn_list_widget.customContextMenuRequested.connect(self.show_connection_menu)
    def show_connection_menu(self, pos):
        item = self.conn_list_widget.itemAt(pos)
        if not item:
            return

        menu = QMenu()
        delete_action = menu.addAction("Bağlantıyı Sil")

        action = menu.exec(self.conn_list_widget.mapToGlobal(pos))
        if action == delete_action:
            text = item.text()  # "a1 <-> a2"
            node1, node2 = [x.strip() for x in text.split("<->")]
            self.remove_connection_logic(node1, node2)

    def _find_all_paths(self, start, end, path=None):
            if path is None: path = []
            path = path + [start]
            if start == end:
                return [path]
            if start not in self.graph:
                return []
            paths = []
            for node in self.graph[start]:
                if node not in path:
                    new_paths = self._find_all_paths(node, end, path)
                    for p in new_paths:
                        if p not in paths:
                            paths.append(p)
            return paths

        # --- ARAYÜZ YÖNETİM FONKSİYONLARI ---
        
    def show_formula_window(self):
            if self.formula_latex:
                if self.formula_window and self.formula_window.isVisible():
                    self.formula_window.activateWindow()
                else:
                    if not self.formula_latex or not any(f.strip() for f in self.formula_latex):
                        QMessageBox.information(self, "Bilgi", "Gösterilecek formül yok.")
                        return

                    self.formula_window = FormulaWindow(self.formula_latex)

                    self.formula_window.show()
            else:
                QMessageBox.information(self, "Bilgi", "Önce 'FORMÜL ÜRET' düğmesine basarak formülleri oluşturmalısınız.")

    
    def on_mode_changed(self, mode_text):

        if mode_text == "Statik Analiz (R)":
            self.analysis_mode = "static"
            self.static_inputs_widget.setVisible(True)
            self.dynamic_inputs_widget.setVisible(False)
            
            self.run_button.setText("FORMÜL ÜRET & HESAPLA")

        elif mode_text == "Dinamik Analiz R(t)":
            self.analysis_mode = "dynamic"
            self.static_inputs_widget.setVisible(False)
            self.dynamic_inputs_widget.setVisible(True)
           
            self.run_button.setText("FORMÜL ÜRET & R(t) ÇİZ")
            self.on_dist_changed()


        else:  # --- Monte Carlo ---
            self.analysis_mode = "montecarlo"
            self.static_inputs_widget.setVisible(False)
            self.dynamic_inputs_widget.setVisible(True)

            self.run_button.setText("Monte Carlo ÇALIŞTIR")
            self.on_dist_changed()

    # --- EKSİK FONKSİYONLARIN SONU ---    
        
    def refresh_left_panel(self):
        # 🔴 HER ŞEYİ TEMİZLE
        self.comp_list_widget.clear()
        self.conn_list_widget.clear()
        self.from_node_selector.clear()
        self.to_node_selector.clear()

        # 🧩 BİLEŞEN LİSTESİ
        for cname, cdata in self.components.items():
            if cdata["dist"] == "static":
                text = f"{cname} (R={cdata['R']:.4f})"
            else:
                params = ", ".join(f"{k}={v}" for k, v in cdata["params"].items())
                text = f"{cname} ({cdata['dist']}: {params})"

            self.comp_list_widget.addItem(text)

        # 🔗 NODE SELECTOR
        all_nodes = ["Start", "End"] + list(self.components.keys()) + sorted(self.junctions)
        self.from_node_selector.addItems(all_nodes)
        self.to_node_selector.addItems(all_nodes)

        # 🔗 BAĞLANTILAR
        for a, b in self.connections:
            self.conn_list_widget.addItem(f"{a} <-> {b}")


    def add_component(self):
        if not hasattr(self, "scene") or self.scene is None:
            QMessageBox.warning(
                self,
                "Aktif Model Yok",
                "Lütfen önce bir model yükleyin veya aktif bir sekme seçin."
            )
            return
        

        comp_name = self.comp_name_input.text().strip()

        if not comp_name or comp_name in self.components or comp_name in ["Start", "End"]:
            QMessageBox.warning(self, "Hata", "Geçersiz veya mevcut bileşen adı.")
            return

        try:
            # --- STATİK ---
            if self.analysis_mode == "static":
                R = float(self.comp_reli_input.text())
                if not (0 <= R <= 1):
                    raise ValueError

                self.components[comp_name] = {
                    "dist": "static",
                    "R": R
                }

                display_text = f"{comp_name} (R={R:.4f})"

            # --- DİNAMİK ---
            else:
                dist = self.dist_selector.currentText()
                params = {}

                dist = self.dist_selector.currentText()
                config = DISTRIBUTIONS[dist]

                for (prm, (_, edit)) in zip(config["params"], self.param_inputs):
                    params[prm["key"]] = float(edit.text())

                self.components[comp_name] = {
                    "dist": dist,
                    "params": params
                }

                param_str = ", ".join(f"{k}={v}" for k, v in params.items())
                display_text = f"{comp_name} ({dist}: {param_str})"


        except Exception:
            QMessageBox.warning(self, "Hata", "Parametreleri doğru giriniz.")
            return

        # Listeye ekle
        self.comp_list_widget.addItem(display_text)
        self.update_connection_selectors()

        # Node çiz
        x = 300 + (len(self.components) % 4) * 120
        y = 100 + (len(self.components) // 4) * 80
        self.draw_node(comp_name, x, y, QColor("darkblue"), is_component=True)

        # Alanları temizle
        self.comp_name_input.clear()
        self.comp_reli_input.clear()

    
    # --- BAĞLANTI VE KAVŞAK MANTIĞI ---

    def add_connection(self):
        if not hasattr(self, "scene") or self.scene is None:
            QMessageBox.warning(
                self,
                "Aktif Model Yok",
                "Lütfen önce bir model yükleyin."
            )
            return

        node1 = self.from_node_selector.currentText()
        node2 = self.to_node_selector.currentText()
        self.add_connection_logic(node1, node2)

    def add_connection_logic(self, node1, node2):
        """ İki düğüm arasına bağlantı ekleyen çekirdek mantık """
        if not node1 or not node2 or node1 == node2:
            QMessageBox.warning(self, "Hata", "Geçersiz bağlantı."); return False
        
        edge_key = tuple(sorted((node1, node2)))
        if edge_key in self.edge_items: 
            QMessageBox.information(self, "Bilgi", "Bu bağlantı zaten mevcut.")
            return False

        if node1 not in self.graph: self.graph[node1] = []
        if node2 not in self.graph: self.graph[node2] = []
        self.graph[node1].append(node2)
        self.graph[node2].append(node1)
        
        display_text = f"{node1} <-> {node2}"
        self.conn_list_widget.addItem(display_text)
        if hasattr(self, "connections"):
            self.connections.append((node1, node2))

        line_item = self.draw_edge(node1, node2)
        if line_item:
            self.edge_items[edge_key] = line_item
        
        return True

    def remove_connection_logic(self, node1, node2):
        """ İki düğüm arasındaki bağlantıyı (mantık ve görsel) siler """
        edge_key = tuple(sorted((node1, node2)))
        
        if edge_key in self.edge_items:
            line_item = self.edge_items.pop(edge_key)
            self.scene.removeItem(line_item)
        
        if node1 in self.graph and node2 in self.graph[node1]:
            self.graph[node1].remove(node2)
        if node2 in self.graph and node1 in self.graph[node2]:
            self.graph[node2].remove(node1)
            
        display_text = f"{node1} <-> {node2}"
        display_text_rev = f"{node2} <-> {node1}"
        for i in range(self.conn_list_widget.count()):
            item = self.conn_list_widget.item(i)
            if item.text() == display_text or item.text() == display_text_rev:
                self.conn_list_widget.takeItem(i)
                break
        if hasattr(self, "connections"):
            if (node1, node2) in self.connections:
                self.connections.remove((node1, node2))
            elif (node2, node1) in self.connections:
                self.connections.remove((node2, node1))
                
    def split_line(self, node1, node2, click_pos):
        """
        'Tele Sağ Tıklama' fonksiyonu.
        Otomatik olarak teli böler ve bir kavşak (junction) ekler.
        """
        j_name = f"j{self.junction_count}"
        self.junction_count += 1
        self.junctions.add(j_name)
        
        self.draw_node(j_name, click_pos.x(), click_pos.y(), QColor("darkgrey"), is_component=False)
        self.update_connection_selectors() # Kavşağı listelere ekle
        
        self.remove_connection_logic(node1, node2)
        
        self.add_connection_logic(node1, j_name)
        self.add_connection_logic(j_name, node2)
    def remove_component(self, comp_name):
        # 1) Bu bileşene bağlı tüm bağlantıları bul
        connections_to_remove = [
            (a, b) for (a, b) in list(self.connections)
            if a == comp_name or b == comp_name
        ]

        # 2) Önce bağlantıları temizle
        for a, b in connections_to_remove:
            self.remove_connection_logic(a, b)

        # 3) Sahnedeki node'u sil
        if comp_name in self.node_items:
            self.scene.removeItem(self.node_items[comp_name])
            del self.node_items[comp_name]

        # 4) Modelden sil
        self.components.pop(comp_name, None)
        self.node_positions.pop(comp_name, None)

        # 5) Sol paneli yenile
        self.refresh_left_panel()

        print(f"[INFO] Component '{comp_name}' removed.")
            
    def update_node_position(self, node_name, new_pos):
        width = 80 if (node_name in self.components or node_name in ['Start', 'End']) else 30
        height = 40 if (node_name in self.components or node_name in ['Start', 'End']) else 30
        
        center_pos = QPointF(new_pos.x() + width / 2, new_pos.y() + height / 2)
        self.node_positions[node_name] = center_pos

        for edge_key, line_item in self.edge_items.items():
            if node_name in edge_key:
                node1, node2 = edge_key
                try:
                    p1_center = self.node_positions[node1]
                    p2_center = self.node_positions[node2]
                    
                    V = p2_center - p1_center; L = math.sqrt(V.x()**2 + V.y()**2)
                    if L == 0: line_item.setLine(0,0,0,0); continue
                    
                    U = V / L
                    clip1 = 42 if (node1 in self.components or node1 in ['Start', 'End']) else 17
                    clip2 = 42 if (node2 in self.components or node2 in ['Start', 'End']) else 17
                    
                    if L <= clip1 + clip2: 
                        line_item.setLine(0,0,0,0); continue
                        
                    p1_clipped = p1_center + U * clip1
                    p2_clipped = p2_center - U * clip2
                    line_item.setLine(p1_clipped.x(), p1_clipped.y(), p2_clipped.x(), p2_clipped.y())
                except Exception as e:
                    pass 

    def reset_model(self):
        self.components = {}; self.graph = {}; self.node_positions = {}
        self.node_items = {}; self.edge_items = {} 
        self.junctions = set(); self.junction_count = 1
        self.scene.clear(); self.comp_list_widget.clear(); self.conn_list_widget.clear()
        self.result_label.setText("Sistem Güvenirliği: - (Model değişti)")
        
        self.formula_latex = None 
        
        self.show_formula_button.setEnabled(True)
        self.initialize_scene()

    # --- Kalan Yardımcı Fonksiyonlar ---
    def create_separator(self, text):
        label = QLabel(text); label.setFont(QFont("Arial", 10, QFont.Weight.Bold)); label.setStyleSheet("margin-top: 8px; margin-bottom: 2px;")
        return label
    def initialize_scene(self):
        scene_w = self.view.width()
        scene_h = self.view.height()

        # Start solda ortalanmış
        start_x = scene_w * 0.15
        start_y = scene_h * 0.5

        # End sağda ortalanmış
        end_x = scene_w * 0.85
        end_y = scene_h * 0.5

        self.draw_node('Start', start_x, start_y, QColor("darkgreen"), is_component=False)
        self.draw_node('End', end_x, end_y, QColor("darkred"), is_component=False)

        self.update_connection_selectors()

        
    def draw_node(self, name, x, y, color, is_component=True):
        if is_component:
             node_item = DraggableNode(0, 0, 80, 40, node_name=name, main_window=self) 
             width, height = 80, 40
        else:
             node_item = DraggableNode(0, 0, 30, 30, node_name=name, main_window=self)
             width, height = 30, 30
             
        node_item.setBrush(QBrush(color)); node_item.setPen(QPen(Qt.GlobalColor.black))
        text_item = QGraphicsTextItem(name, node_item); text_item.setDefaultTextColor(QColor("white")); text_item.setFont(QFont("Arial", 9, QFont.Weight.Bold))
        text_rect = text_item.boundingRect(); text_item.setPos((width - text_rect.width()) / 2, (height - text_rect.height()) / 2); self.scene.addItem(node_item)
        
        node_item.setPos(x - width/2, y - height/2) 
        
        self.node_positions[name] = QPointF(x, y) # Pozisyonu MERKEZ olarak kaydet
        self.node_items[name] = node_item; 
        
        self.update_node_position(name, node_item.pos())
        
        return node_item
        
    def draw_edge(self, from_node, to_node):
        """ Düz çizgi (ClickableLineItem) kullanır """
        try:
            pos1 = self.node_positions[from_node]; pos2 = self.node_positions[to_node]
            
            p1_center = pos1; p2_center = pos2
            
            V = p2_center - p1_center; L = math.sqrt(V.x()**2 + V.y()**2)
            if L == 0: return None
            
            U = V / L
            clip1 = 42 if (from_node in self.components or from_node in ['Start', 'End']) else 17
            clip2 = 42 if (to_node in self.components or to_node in ['Start', 'End']) else 17
            
            if L <= clip1 + clip2: return None 
            
            p1_clipped = p1_center + U * clip1
            p2_clipped = p2_center - U * clip2
            
            line = ClickableLineItem(from_node, to_node, self, 
                                     p1_clipped.x(), p1_clipped.y(), 
                                     p2_clipped.x(), p2_clipped.y())
            
            self.scene.addItem(line)
            return line
        except KeyError: 
            print(f"Hata: {from_node} veya {to_node} için pozisyon bulunamadı.")
            return None
            
    def update_connection_selectors(self):
        self.from_node_selector.clear(); self.to_node_selector.clear()
        all_nodes = ['Start', 'End'] + list(self.components.keys()) + sorted(list(self.junctions))
        self.from_node_selector.addItems(all_nodes)
        self.to_node_selector.addItems(all_nodes)
    

    def save_scenario_to_json(self, scenario_name, params, t_values, system_rt, comp_rt, T_sys_samples=None):

        """
        Tek bir analizi JSON dosyasına ekler.
        scenario_name : kullanıcının verdiği isim (örn: 'test1')
        params        : bileşen parametreleri (self.components)
        t_values      : numpy array (zaman noktaları)
        system_rt     : numpy array (Sistem R(t))
        comp_rt       : dict {comp_name: numpy array}
        """

        filename = "scenarios.json"

        # JSON kayıt formatı
        entry = {
    "name": scenario_name,
    "analysis_mode": self.analysis_mode,
    "params": params,
    "t_values": list(map(float, t_values)),
    "system_rt": list(map(float, system_rt)),
    "components_rt": {k: list(map(float, v)) for k, v in comp_rt.items()}
}

        if T_sys_samples is not None:
            entry["T_sys_samples"] = T_sys_samples

        # Eski dosyayı oku
        if os.path.exists(filename):
            with open(filename, "r") as f:
                try:
                    existing = json.load(f)
                except json.JSONDecodeError:
                    existing = {"scenarios": []}
        else:
            existing = {"scenarios": []}

        # Aynı isimde senaryo varsa → üzerine yaz
        existing["scenarios"] = [
            sc for sc in existing.get("scenarios", [])
            if sc.get("name") != scenario_name
        ]

        existing["scenarios"].append(entry)

        # Yaz
        with open(filename, "w") as f:
            json.dump(existing, f, indent=4)

        print(f"[JSON] Senaryo kaydedildi: {scenario_name}")


    def plot_from_json(self, scenario_names):
        filename = "scenarios.json"

        if not os.path.exists(filename):
            QMessageBox.warning(self, "Hata", "Önce senaryo kaydedilmesi gerekiyor!")
            return

        with open(filename, "r") as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                QMessageBox.warning(self, "Hata", "scenarios.json okunamadı.")
                return

        t_values = None
        plot_data = {}

        for sc in data.get("scenarios", []):
            name = sc.get("name")
            if name in scenario_names:

                # Eksik veri varsa atla
                if "system_rt" not in sc or "t_values" not in sc:
                    print(f"[UYARI] Senaryo bozuk atlanıyor: {name}")
                    continue

                if t_values is None:
                    t_values = np.array(sc["t_values"], dtype=float)

                plot_data[name] = np.array(sc["system_rt"], dtype=float)

        if not plot_data:
            QMessageBox.warning(self, "Hata", "Seçilen senaryolar bulunamadı veya eksik.")
            return

        # Grafik çiz
        self.plot_window = PlotWindow(t_values, plot_data)
        self.plot_window.show()
    

    
    def run_main_button(self):
        if self.analysis_mode == "static":
            self.run_analysis()  # mevcut analitik statik fonksiyonun

        elif self.analysis_mode == "dynamic":
            self.run_analysis()  # mevcut analitik dinamik fonksiyonun

        elif self.analysis_mode == "montecarlo":
            self.run_monte_carlo_gui()


    def clean_latex(self, text):
    
    
        return text
    def save_model(self, filename):
        """
        Sadece bileşenleri, kavşakları, bağlantıları ve
        Start/End HARİÇ düğüm pozisyonlarını kaydeder.
        """
        data = {
            "analysis_mode": getattr(self, "analysis_mode", "static"),
            "components": self.components,
            "junctions": list(self.junctions),
            "node_positions": {
                name: [
                    float(self.node_positions[name].x()),
                    float(self.node_positions[name].y())
                ]
                for name in self.node_positions
                if name not in ["Start", "End"]        # <<< ÖNEMLİ
            },
            "connections": [list(edge) for edge in self.edge_items.keys()]
        }

        with open(filename, "w") as f:
            json.dump(data, f, indent=4)

        print("[MODEL] Kaydedildi:", filename)

    def load_model(self, filename):
        # 1️⃣ Aktif sekmenin state’ini al
        tab = self.tab_widget.currentWidget()
        state = tab.model_state

        # 2️⃣ JSON oku
        with open(filename, "r") as f:
            data = json.load(f)

        # 3️⃣ STATE’İ TEMİZLE
        state["components"].clear()
        state["junctions"].clear()
        state["connections"].clear()
        state["graph"].clear()
        state["node_positions"].clear()
        state["node_items"].clear()
        state["edge_items"].clear()

        # 4️⃣ STATE’E JSON VERİSİNİ YAZ
        state["components"].update(copy.deepcopy(data.get("components", {})))
        state["junctions"].update(set(data.get("junctions", [])))
        state["node_positions"].update(data.get("node_positions", {}))

        # 5️⃣ SCENE TEMİZLE + START / END
        self.scene.clear()
        self.initialize_scene()

        # 6️⃣ BİLEŞENLERİ ÇİZ
        for cname, cdata in state["components"].items():
            if cname in ["Start", "End"]:
                continue

            if cname in state["node_positions"]:
                x, y = state["node_positions"][cname]
            else:
                x, y = 300, 200

            self.draw_node(cname, x, y, QColor("darkblue"), is_component=True)

        # 7️⃣ KAVŞAKLARI ÇİZ
        for jname in state["junctions"]:
            if jname in state["node_positions"]:
                x, y = state["node_positions"][jname]
                self.draw_node(jname, x, y, QColor("darkgrey"), is_component=False)

        # 8️⃣ BAĞLANTILARI KUR
        for n1, n2 in data.get("connections", []):
            self.add_connection_logic(n1, n2)

        # 9️⃣ SOL PANEL + STATE SENKRON
        self.on_tab_changed(self.tab_widget.currentIndex())

        # 10️⃣ UI DURUM
        self.result_label.setText("Sistem Güvenirliği: - (Model yüklendi)")
        print("[MODEL] Yüklendi:", filename)

    def save_model_dialog(self):
        # Kaydetme penceresi aç
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Modeli Kaydet",
            "",
            "JSON Files (*.json)"
        )

        if not filename:
            return

        # uzantı yoksa ekle
        if not filename.lower().endswith(".json"):
            filename += ".json"

        try:
            self.save_model(filename)
            QMessageBox.information(self, "Başarılı", f"Model kaydedildi:\n{filename}")
        except Exception as e:
            QMessageBox.critical(self, "Kaydetme Hatası", f"Model kaydedilemedi:\n{e}")

    def load_model_dialog(self):
        
        filenames, _ = QFileDialog.getOpenFileNames(
            self,
            "Bir veya Daha Fazla Model JSON Seç",
            "",
            "JSON Files (*.json)"
        )

        if not filenames:
            return

        # === ÇOKLU MODEL ANALİZİ ===
        self.loaded_models = {}

        for filepath in filenames:
            model_name = os.path.basename(filepath)
            print(f"[BATCH] Model yükleniyor: {model_name}")

            # === YENİ SEKME OLUŞTUR ===
            tab = QWidget()
            tab_layout = QVBoxLayout(tab)

            scene = QGraphicsScene()
            scene.setBackgroundBrush(QBrush(QColor(230, 230, 230)))

            view = QGraphicsView(scene)
            view.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
            view.setRenderHint(QPainter.RenderHint.Antialiasing)

            tab_layout.addWidget(view)
            tab.setLayout(tab_layout)

            self.tab_widget.addTab(tab, model_name)
            tab.model_state = {
    "components": {},
    "junctions": set(),
    "connections": [],
    "graph": {},
    "node_positions": {},
    "node_items": {},
    "edge_items": {}
}
            self.tab_widget.setCurrentWidget(tab)
            self.on_tab_changed(self.tab_widget.currentIndex())


            # 👉 AKTİF SAHNEYİ BU SEKME YAP
            self.scene = scene
            self.view = view
            # 1) Modeli yükle
            self.load_model(filepath)


            # 3) Sonucu sakla
            if hasattr(self, "last_results"):
                self.loaded_models[model_name] = {
                    "t": self.last_results["CurrentSystem"]["t"],
                    "R": self.last_results["CurrentSystem"]["R"]
                }

        # === KRİTİK BÖLGE KARŞILAŞTIRMA RAPORLARI ===
        if len(self.loaded_models) >= 2:
            intervals = plot_critical_intervals(self.loaded_models)
            plot_critical_slopes(intervals)
            

            df = critical_summary_table(intervals)
            print("\n=== KRİTİK ANALİZ ÖZET TABLOSU ===")
            print(df)
        else:
            QMessageBox.information(
                self,
                "Bilgi",
                "Kritik karşılaştırma için en az 2 model yüklenmelidir."
            )

    from monte_carlo import run_monte_carlo
    def _get_ccf_config(self):
        """
        CCF (Common Cause Failure) ayarlarını döndürür.
        Aktif değilse None döner.
        """
        if not self.ccf_checkbox.isChecked():
            return None

        try:
            beta = float(self.ccf_beta_input.text())
            if not (0.0 <= beta <= 1.0):
                raise ValueError

            # Sadece Exponential bileşenlerden lambda al
            lambdas = [
                d["params"]["lambda"]
                for d in self.components.values()
                if d["dist"] == "Exponential" and "lambda" in d["params"]
            ]

            if not lambdas:
                return None

            return (beta, lambdas)

        except Exception:
            QMessageBox.warning(
                self,
                "CCF Hatası",
                "CCF β değeri 0–1 arasında olmalı ve en az bir Exponential bileşen olmalı."
            )
            return None

    def run_monte_carlo_gui(self):
        tab = self.tab_widget.currentWidget()
        if "component_paths" not in tab.model_state:
            QMessageBox.warning(
                self,
                "Monte Carlo Uyarısı",
                "Monte Carlo simülasyonu için önce Dinamik Analiz (R(t)) çalıştırılmalıdır.\n\n"
                "Sebep: Monte Carlo, analitik olarak elde edilen minimal yol kümelerini kullanır."
            )
            return
        self.component_paths = tab.model_state.get("component_paths", [])
        T_sys, t_vals, R_mc = run_monte_carlo(
        components=self.components,
        component_paths=self.component_paths,
        N=self.mc_spinbox.value(),
        t_max=self.t_max_input.value(),
        ccf=self._get_ccf_config()
    )

        # Histogram
        self.mc_hist_window = HistogramWindow(T_sys)
        self.mc_hist_window.setWindowTitle("Monte Carlo Sistem Ömrü Histogramı")
        self.mc_hist_window.show()

        # 🔥 R(t) grafiği
        plot_data = {
            "Monte Carlo R(t)": R_mc
        }
        self.plot_window = PlotWindow(t_vals, plot_data)
        self.plot_window.show()

    def resizeEvent(self, event):
        super().resizeEvent(event)

        # Start & End zaten varsa yeniden konumlandır
        if "Start" in self.node_items and "End" in self.node_items:
            scene_w = self.view.width()
            scene_h = self.view.height()

            start_x = scene_w * 0.15
            start_y = scene_h * 0.5
            end_x   = scene_w * 0.85
            end_y   = scene_h * 0.5

            self.node_items["Start"].setPos(start_x - 15, start_y - 15)
            self.node_positions["Start"] = QPointF(start_x, start_y)

            self.node_items["End"].setPos(end_x - 15, end_y - 15)
            self.node_positions["End"] = QPointF(end_x, end_y)

            # bağlantıları güncelle
            self.update_node_position("Start", self.node_items["Start"].pos())
            self.update_node_position("End", self.node_items["End"].pos())
    def on_component_selected(self, row):
        if row >= 0:
            self.update_comp_button.setEnabled(True)

            # Seçili bileşenin adını inputlara yaz
            item_text = self.comp_list_widget.item(row).text()
            comp_name = item_text.split()[0]  # a1, a2 vs

            self.comp_name_input.setText(comp_name)

            data = self.components.get(comp_name)
            if not data:
                return

            if data["dist"] == "static":
                self.comp_reli_input.setText(str(data["R"]))
            else:
                self.dist_selector.setCurrentText(data["dist"])
                self.on_dist_changed()
                for (label, edit), (_, val) in zip(
                    self.param_inputs, data["params"].items()
                ):
                    edit.setText(str(val))
        else:
            self.update_comp_button.setEnabled(False)
    def update_selected_component(self):
        row = self.comp_list_widget.currentRow()
        if row < 0:
            return

        item = self.comp_list_widget.item(row)
        comp_name = item.text().split()[0]

        try:
            if self.analysis_mode == "static":
                R = float(self.comp_reli_input.text())
                if not (0 <= R <= 1):
                    raise ValueError

                self.components[comp_name] = {
                    "dist": "static",
                    "R": R
                }
                item.setText(f"{comp_name} (R={R:.4f})")

            else:
                dist = self.dist_selector.currentText()
                config = DISTRIBUTIONS[dist]
                params = {}

                for (prm, (_, edit)) in zip(config["params"], self.param_inputs):
                    params[prm["key"]] = float(edit.text())

                self.components[comp_name] = {
                    "dist": dist,
                    "params": params
                }

                param_str = ", ".join(f"{k}={v}" for k, v in params.items())
                item.setText(f"{comp_name} ({dist}: {param_str})")

            QMessageBox.information(self, "Başarılı", f"{comp_name} güncellendi.")

        except Exception:
            QMessageBox.warning(self, "Hata", "Parametreleri doğru giriniz.")
    def _compute_mttf_for_components(self, components_backup):
        """
        Verilen bileşen konfigürasyonu için
        SADECE MTTF hesaplar (grafik çizmez)
        """
        # === minimal yolları bul ===
        paths_nodes = self._find_all_paths("Start", "End")
        component_paths = [
            frozenset([n for n in p if n in components_backup])
            for p in paths_nodes
        ]
        component_paths = list(dict.fromkeys(component_paths))

        t_max = self.t_max_input.value()
        t_safe = np.linspace(1e-6, t_max, 400)

        t = sympy.symbols("t", positive=True)
        modules = ["numpy", {"exp": np.exp, "log": np.log, "sqrt": np.sqrt, "erfc": scipy_erfc}]

        # === path R(t) ===
        path_rts = []

        for pset in component_paths:
            rt_path = np.ones_like(t_safe)

            for cname in pset:
                data = components_backup[cname]
                if data["dist"] == "static":
                    rt_c = np.ones_like(t_safe) * data["R"]
                else:
                    conf = DISTRIBUTIONS[data["dist"]]
                    if conf.get("R_sym") is not None:
                        f = sympy.lambdify(t, conf["R_sym"](t, data["params"]), modules)
                        rt_c = f(t_safe)
                    else:
                        rt_c = conf["R_num"](t_safe, data["params"])


                if conf.get("R_sym") is not None:
                    f = sympy.lambdify(t, conf["R_sym"](t, data["params"]), modules)
                    rt_c = f(t_safe)
                else:
                    rt_c = conf["R_num"](t_safe, data["params"])

                rt_c = np.nan_to_num(rt_c, 0.0, 0.0, 0.0)
                rt_path *= rt_c

            path_rts.append(rt_path)

        # === inclusion–exclusion ===
        system_r = np.zeros_like(t_safe)
        for k in range(1, len(path_rts) + 1):
            for comb in combinations(range(len(path_rts)), k):
                prod = np.ones_like(t_safe)
                for i in comb:
                    prod *= path_rts[i]
                system_r += ((-1) ** (k + 1)) * prod

        system_r = np.clip(system_r, 0, 1)
        mttf = np.trapz(system_r, t_safe)

        return float(mttf)

    def run_sensitivity_analysis(self):
        """
        Her bileşenin parametresini %±10 değiştirir,
        MTTF farkını hesaplar (ΔMTTF)
        """
        base_components = copy.deepcopy(self.components)
        mttf_base = self._compute_mttf_for_components(base_components)

        sensitivity_results = {}

        for cname, cdata in base_components.items():
            if cdata["dist"] == "static":
                continue  # statik için anlamlı değil

            params = cdata["params"]

            # sadece ilk parametreyi oynatıyoruz (literatürde yaygın)
            key = list(params.keys())[0]
            original = params[key]

            if original <= 0:
                continue

            # +10%
            params[key] = original * 1.10
            mttf_plus = self._compute_mttf_for_components(base_components)

            # -10%
            params[key] = original * 0.90
            mttf_minus = self._compute_mttf_for_components(base_components)

            # geri al
            params[key] = original

            # en büyük etkiyi al
            delta = max(
                abs(mttf_plus - mttf_base),
                abs(mttf_minus - mttf_base)
            )

            sensitivity_results[cname] = delta

        return mttf_base, sensitivity_results


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())