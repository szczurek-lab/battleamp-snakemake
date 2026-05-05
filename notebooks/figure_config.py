"""Shared configuration for all BATTLE-AMP manuscript figures."""
from pathlib import Path
import os
import sys
import json
import logging
import string
import csv
import numpy as np
import pandas as pd
import warnings

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("battleamp")

# Save control: 'pdf', 'png', 'both', or None
SAVE = "both"

# Paths (scripts live in battleamp-snakemake/notebooks/)
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
FIGURE_DIR   = PROJECT_ROOT / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

DATA_ROOT    = PROJECT_ROOT / "data"
TASKS_DIR    = PROJECT_ROOT / "tasks"
RESULTS_ROOT = PROJECT_ROOT / "results"
CLF_FILE     = RESULTS_ROOT / "aggregated" / "classification_results.tsv"
REG_FILE     = RESULTS_ROOT / "aggregated" / "regression_results.tsv"
INFERENCE_DIR = RESULTS_ROOT / "inference"
MD_DIR       = PROJECT_ROOT / "md"
CACHE_DIR    = PROJECT_ROOT / "figure4_cache"

# Journal constraints
MAX_WIDTH  = 6.5
MAX_HEIGHT = 8.0
DPI        = 600

# Font sizes (all >= 6 pt journal minimum)
FS_TICK    = 6
FS_LABEL   = 7
FS_PANEL   = 10
FS_ANNOT   = 6
FS_TITLE   = 7
FS_LEGEND  = 6

# Standardised colour palette
ACTIVE_COLOR   = "#F7CF8B"
INACTIVE_COLOR = "#66BDBA"

RANDOM_COLOR    = "#7E7F9A"
SHUFFLED_COLOR  = "#FF6542"
REALISTIC_COLOR = "#721817"

AMP_COLOR = "#8B5E3C"
GA_COLOR  = "#2B6A99"

ZP_NEG_COLOR = "#5B8DB8"
ZP_POS_COLOR = "#B8725B"

# Model-type colours (Okabe-Ito, colorblind-friendly)
CLF_CLR       = "#E69F00"
ACT_CLR       = "#CC79A7"
REG_CLR       = "#0072B2"
CLF_CLR_LIGHT = "#F3CF80"
ACT_CLR_LIGHT = "#E4BCCF"
REG_CLR_LIGHT = "#80B9D9"

DATASET_PAL = {
    "AMP":             AMP_COLOR,
    "GeneralActivity": GA_COLOR,
    "Random":          RANDOM_COLOR,
    "Shuffled":        SHUFFLED_COLOR,
    "Realistic":       REALISTIC_COLOR,
}

BG         = "#F0F0F0"
GRID_CLR   = "#cccccc"
TEXT_CLR   = "#000000"
EXCLUDED_ZONE = "#E0E0E0"

# ── matplotlib ────────────────────────────────────────────────────────────────
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns

matplotlib.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.default":   "regular",
    "font.size":          FS_LABEL,
    "axes.labelsize":     FS_LABEL,
    "axes.titlesize":     FS_TITLE,
    "xtick.labelsize":    FS_TICK,
    "ytick.labelsize":    FS_TICK,
    "legend.fontsize":    FS_LEGEND,
    "axes.linewidth":     0.4,
    "xtick.major.width":  0.4,
    "ytick.major.width":  0.4,
    "xtick.major.size":   2,
    "ytick.major.size":   2,
    "axes.spines.right":  False,
    "axes.spines.top":    False,
    "figure.facecolor":   "white",
    "axes.facecolor":     BG,
    "axes.grid":          False,
    "text.color":         TEXT_CLR,
    "axes.labelcolor":    TEXT_CLR,
    "xtick.color":        TEXT_CLR,
    "ytick.color":        TEXT_CLR,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
})


def style_ax(ax):
    ax.set_facecolor(BG)
    ax.yaxis.grid(True, color=GRID_CLR, linewidth=0.3, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig, stem):
    if SAVE is None:
        return
    fmt = SAVE.lower() if isinstance(SAVE, str) else None
    kw = dict(bbox_inches="tight", pad_inches=0.02, dpi=DPI, facecolor="white")
    if fmt in ("pdf", "both"):
        p = FIGURE_DIR / f"{stem}.pdf"
        fig.savefig(str(p), **kw)
        print(f"Saved {p}")
    if fmt in ("png", "both"):
        p = FIGURE_DIR / f"{stem}.png"
        fig.savefig(str(p), **kw)
        print(f"Saved {p}")

# ── model nomenclature ────────────────────────────────────────────────────────
from collections import OrderedDict

CLASSIFIERS = [
    "ampeppy", "amplify", "ampredmfa", "ampscanner",
    "hydramp-amp-classifier", "sensexamp-classifier",
]
ACTIVITY_AWARE = ["hydramp-mic-classifier"]
REGRESSORS = [
    "ampredictor",
    "apex-abaumannii", "apex-ecoli", "apex-kpneumoniae",
    "apex-min", "apex-paeruginosa", "apex-saureus",
    "deep-amp-cnn-gramneg", "deep-amp-cnn-grampos",
    "deep-amp-lstm-gramneg", "deep-amp-lstm-grampos",
    "mbc-attention",
    "sensexamp-ecoli", "sensexamp-saureus",
]
ALL_MODELS = CLASSIFIERS + ACTIVITY_AWARE + REGRESSORS

MODEL_DISPLAY = {
    "hydramp-amp-classifier":  r"HydrAMP$_{\mathrm{AMP}}$",
    "ampscanner":              r"AMPScanner$_{\mathrm{v2}}$",
    "amplify":                 "AMPlify",
    "sensexamp-classifier":    r"SenseXAMP$_{\mathrm{clf}}$",
    "ampeppy":                 "amPEPpy",
    "ampredmfa":               "AMPpred-MFA",
    "hydramp-mic-classifier":  r"HydrAMP$_{\mathrm{Ec}}$",
    "mbc-attention":           "MBC-Attention",
    "ampredictor":             "AMPredictor",
    "sensexamp-ecoli":         r"SenseXAMP$_{\mathrm{Ec}}$",
    "sensexamp-saureus":       r"SenseXAMP$_{\mathrm{Sa}}$",
    "deep-amp-lstm-gramneg":   r"DeepAMP$_{\mathrm{LSTM,\,G-}}$",
    "deep-amp-lstm-grampos":   r"DeepAMP$_{\mathrm{LSTM,\,G+}}$",
    "deep-amp-cnn-gramneg":    r"DeepAMP$_{\mathrm{CNN,\,G-}}$",
    "deep-amp-cnn-grampos":    r"DeepAMP$_{\mathrm{CNN,\,G+}}$",
    "apex-ecoli":              r"APEX$_{\mathrm{Ec}}$",
    "apex-saureus":            r"APEX$_{\mathrm{Sa}}$",
    "apex-min":                r"APEX$_{\mathrm{min}}$",
    "apex-abaumannii":         r"APEX$_{\mathrm{Ab}}$",
    "apex-paeruginosa":        r"APEX$_{\mathrm{Pa}}$",
    "apex-kpneumoniae":        r"APEX$_{\mathrm{Kp}}$",
}

CLF_FIG4 = OrderedDict([
    ("ampeppy",                 "amPEPpy"),
    ("amplify",                 "AMPlify"),
    ("ampredmfa",               "AMPpred-MFA"),
    ("ampscanner",              r"AMPScanner$_{\mathrm{v2}}$"),
    ("hydramp-amp-classifier",  r"HydrAMP$_{\mathrm{AMP}}$"),
    ("sensexamp-classifier",    r"SenseXAMP$_{\mathrm{clf}}$"),
    ("hydramp-mic-classifier",  r"HydrAMP$_{\mathrm{Ec}}$"),
])

REG_FIG4 = OrderedDict([
    ("ampredictor",             "AMPredictor"),
    ("apex-abaumannii",         r"APEX$_{\mathrm{Ab}}$"),
    ("apex-ecoli",              r"APEX$_{\mathrm{Ec}}$"),
    ("apex-kpneumoniae",        r"APEX$_{\mathrm{Kp}}$"),
    ("apex-min",                r"APEX$_{\mathrm{min}}$"),
    ("apex-paeruginosa",        r"APEX$_{\mathrm{Pa}}$"),
    ("apex-saureus",            r"APEX$_{\mathrm{Sa}}$"),
    ("deep-amp-cnn-gramneg",    r"DeepAMP$_{\mathrm{CNN,\,G-}}$"),
    ("deep-amp-cnn-grampos",    r"DeepAMP$_{\mathrm{CNN,\,G+}}$"),
    ("deep-amp-lstm-gramneg",   r"DeepAMP$_{\mathrm{LSTM,\,G-}}$"),
    ("deep-amp-lstm-grampos",   r"DeepAMP$_{\mathrm{LSTM,\,G+}}$"),
    ("mbc-attention",           "MBC-Attention"),
    ("sensexamp-ecoli",         r"SenseXAMP$_{\mathrm{Ec}}$"),
    ("sensexamp-saureus",       r"SenseXAMP$_{\mathrm{Sa}}$"),
])

HEATMAP_ORDER = list(CLF_FIG4.keys()) + list(REG_FIG4.keys())

SCATTER_DISPLAY = {
    "mbc-attention":         "MBC-Att.",
    "ampredictor":           "AMPred.",
    "sensexamp-ecoli":       r"SX$_{\rm Ec}$",
    "sensexamp-saureus":     r"SX$_{\rm Sa}$",
    "deep-amp-lstm-gramneg": r"DA$_{\rm L,G-}$",
    "deep-amp-lstm-grampos": r"DA$_{\rm L,G+}$",
    "deep-amp-cnn-gramneg":  r"DA$_{\rm C,G-}$",
    "deep-amp-cnn-grampos":  r"DA$_{\rm C,G+}$",
    "apex-ecoli":            r"AP$_{\rm Ec}$",
    "apex-saureus":          r"AP$_{\rm Sa}$",
    "apex-min":              r"AP$_{\rm min}$",
    "apex-abaumannii":       r"AP$_{\rm Ab}$",
    "apex-paeruginosa":      r"AP$_{\rm Pa}$",
    "apex-kpneumoniae":      r"AP$_{\rm Kp}$",
}

ML_MODELS_FIG5 = [
    ("ampredictor",       "ampredictor_MIC",        "ampredictor",         "AMPredictor"),
    ("apex-min",          "apex-min_MIC",           "apex-min",            r"APEX$_{\mathrm{min}}$"),
    ("deep-amp-gram-",    "deep-amp-gram-_MIC",     "deep-amp-lstm-gramneg", r"DeepAMP$_{\mathrm{G-}}$"),
    ("deep-amp-gram+",    "deep-amp-gram+_MIC",     "deep-amp-lstm-grampos", r"DeepAMP$_{\mathrm{G+}}$"),
    ("mbc-attention",     "mbc-attention_MIC",      "mbc-attention",       "MBC-Att."),
    ("sensexamp-ecoli",   "sensexamp-ecoli_MIC",    "sensexamp-ecoli",     r"SenseXAMP$_{\mathrm{Ec}}$"),
    ("sensexamp-saureus", "sensexamp-saureus_MIC",  "sensexamp-saureus",   r"SenseXAMP$_{\mathrm{Sa}}$"),
]


def short_name(m):
    return MODEL_DISPLAY.get(m, m)

def scatter_name(m):
    return SCATTER_DISPLAY.get(m, short_name(m))

def model_type(m):
    if m in CLASSIFIERS:    return "clf"
    if m in ACTIVITY_AWARE: return "act"
    return "reg"

def model_color(m):
    return {"clf": CLF_CLR, "act": ACT_CLR}.get(model_type(m), REG_CLR)

def type_bar_dark(t):
    return {"clf": CLF_CLR, "act": ACT_CLR}.get(t, "#999999")

def type_bar_light(t):
    return {"clf": CLF_CLR_LIGHT, "act": ACT_CLR_LIGHT}.get(t, REG_CLR_LIGHT)


# ── shared data-loading utility ───────────────────────────────────────────────
def _normalise_seq_col(df):
    if "Sequence" in df.columns and "sequence" not in df.columns:
        df = df.rename(columns={"Sequence": "sequence"})
    if "sequence" not in df.columns:
        for c in df.columns:
            if c.lower() == "sequence": df = df.rename(columns={c: "sequence"}); break
    return df
