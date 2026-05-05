"""Figure 2: Activity-based evaluation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

df2 = pd.read_csv(CLF_FILE, sep="\t")
df2 = df2[df2["variant"] != "example-model"].copy()

amp     = df2[df2["task"] == "amp"].set_index("variant")
ga      = df2[df2["task"] == "broad_activity"].set_index("variant")
slay_df = df2[df2["task"] == "slay"].set_index("variant")

ga_models = [m for m in ALL_MODELS if m in ga.index]
N2 = len(ga_models)

def _safe(fr, m, col, d=np.nan):
    return fr.loc[m, col] if m in fr.index else d

def _lr(fr, m):
    tpr = _safe(fr, m, "tpr", 0); fpr = _safe(fr, m, "fpr", 0)
    return tpr / fpr if fpr > 0 else 0.0

types2    = [model_type(m) for m in ga_models]
names2    = [short_name(m) for m in ga_models]
amp_mcc   = np.array([_safe(amp, m, "mcc",          0) for m in ga_models])
amp_fpr   = np.array([_safe(amp, m, "fpr",          0) for m in ga_models])
ga_mcc    = np.array([_safe(ga,  m, "mcc",          0) for m in ga_models])
ga_fpr    = np.array([_safe(ga,  m, "fpr",          0) for m in ga_models])
ga_p100   = np.array([_safe(ga,  m, "precision_at_k",0) for m in ga_models])
ga_lr     = np.array([_lr(ga,  m) for m in ga_models])
slay_p100 = np.array([_safe(slay_df, m, "precision_at_k", 0) for m in ga_models])
slay_lr   = np.array([_lr(slay_df, m) for m in ga_models])
slay_base = slay_df.iloc[0]["n_positive"] / (
    slay_df.iloc[0]["n_positive"] + slay_df.iloc[0]["n_negative"])

colors2 = [model_color(m) for m in ga_models]

COMPACT_DISPLAY = {
    "hydramp-amp-classifier":  r"HydrAMP$_{\rm AMP}$",
    "ampscanner":              r"AMPScanner$_{\rm v2}$",
    "amplify":                 "AMPlify",
    "sensexamp-classifier":    r"SenseXAMP$_{\rm clf}$",
    "ampeppy":                 "amPEPpy",
    "ampredmfa":               "AMPpred-MFA",
    "hydramp-mic-classifier":  r"HydrAMP$_{\rm Ec}$",
    "mbc-attention":           "MBC-Attention",
    "ampredictor":             "AMPredictor",
    "sensexamp-ecoli":         r"SenseXAMP$_{\rm Ec}$",
    "sensexamp-saureus":       r"SenseXAMP$_{\rm Sa}$",
    "deep-amp-lstm-gramneg":   r"DeepAMP$_{\rm L,G-}$",
    "deep-amp-lstm-grampos":   r"DeepAMP$_{\rm L,G+}$",
    "deep-amp-cnn-gramneg":    r"DeepAMP$_{\rm C,G-}$",
    "deep-amp-cnn-grampos":    r"DeepAMP$_{\rm C,G+}$",
    "apex-ecoli":              r"APEX$_{\rm Ec}$",
    "apex-saureus":            r"APEX$_{\rm Sa}$",
    "apex-min":                r"APEX$_{\rm min}$",
    "apex-abaumannii":         r"APEX$_{\rm Ab}$",
    "apex-paeruginosa":        r"APEX$_{\rm Pa}$",
    "apex-kpneumoniae":        r"APEX$_{\rm Kp}$",
}
names2c = [COMPACT_DISPLAY.get(m, m) for m in ga_models]


def _group_rank(m):
    """Primary sort key: 0=classifier, 1=HydrAMP_Ec (activity-aware), 2=regressor."""
    return {"clf": 0, "act": 1, "reg": 2}[model_type(m)]


def _group_separators(order):
    """Y-positions between any two consecutive model groups."""
    seps = []
    prev = _group_rank(ga_models[order[0]])
    for rank, i in enumerate(order[1:], 1):
        g = _group_rank(ga_models[i])
        if g != prev:
            seps.append(rank - 0.5)
        prev = g
    return seps


def _ylabels2(ax, order):
    """Y-tick labels coloured by model type: red=classifier, purple=HydrAMP_Ec, grey=regressor."""
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([names2c[i] for i in order], fontsize=FS_TICK)
    for tick, idx in zip(ax.get_yticklabels(), order):
        tick.set_color(model_color(ga_models[idx]))
    ax.set_ylim(-0.6, len(order) - 0.4)


def dumbbell(ax, vf, vo, order, xlabel, xlim=None):
    for rank, i in enumerate(order):
        c = colors2[i]
        ax.plot([vf[i], vo[i]], [rank, rank], color=GRID_CLR, lw=0.8, zorder=1)
        ax.scatter(vf[i], rank, c=c, s=22, zorder=3, edgecolor="white", lw=0.3)
        ax.scatter(vo[i], rank, facecolors="none", edgecolors=c, s=22, lw=1.0, zorder=3)
    ax.axvline(np.nanmedian(vf), color="#444", lw=1.2, ls="-",       alpha=0.7, zorder=0)
    ax.axvline(np.nanmedian(vo), color="#444", lw=1.2, ls=(0,(5,2)), alpha=0.7, zorder=0)
    for sep_y in _group_separators(order):
        ax.axhspan(sep_y - 0.18, sep_y + 0.18, color="white", zorder=2, lw=0)
    _ylabels2(ax, order)
    ax.set_xlabel(xlabel, fontsize=FS_LABEL)
    if xlim: ax.set_xlim(xlim)
    style_ax(ax)
    ax.xaxis.grid(True, color=GRID_CLR, linewidth=0.3)


def grouped_barh(ax, vd, vl, order, xlabel, xlim=None, log=False, ref=None):
    """solid bars = GeneralActivity, tiny-dot bars = SLAY, both colored by model type."""
    h = 0.32
    for rank, i in enumerate(order):
        c = colors2[i]
        ax.barh(rank + h/2, vl[i], height=h, color=c,
                edgecolor=GRID_CLR, lw=0.3)
        ax.barh(rank - h/2, vd[i], height=h, color=c,
                edgecolor="white", lw=0.3, hatch="//")
    for sep_y in _group_separators(order):
        ax.axhspan(sep_y - 0.22, sep_y + 0.22, color="white", zorder=2, lw=0)
    _ylabels2(ax, order)
    ax.tick_params(axis="y", pad=1)
    ax.set_xlabel(xlabel, fontsize=FS_LABEL)
    if log: ax.set_xscale("log")
    if xlim: ax.set_xlim(xlim)
    if ref is not None: ax.axvline(ref, color="#aaa", lw=0.5, ls="--", zorder=0)
    style_ax(ax)
    ax.xaxis.grid(True, color=GRID_CLR, linewidth=0.3)
    valid_l = [vl[i] for i in order if not np.isnan(vl[i]) and vl[i] > 0]
    valid_d = [vd[i] for i in order if not np.isnan(vd[i]) and vd[i] > 0]
    if valid_l: ax.axvline(np.median(valid_l), color="#444", lw=1.2, ls="-",       alpha=0.7, zorder=0)
    if valid_d: ax.axvline(np.median(valid_d), color="#444", lw=1.2, ls=(0,(1,1)), alpha=0.7, zorder=0)


_alpha = {m: k for k, m in enumerate(ALL_MODELS)}
order_ab = sorted(range(N2), key=lambda i: (-_group_rank(ga_models[i]), -_alpha.get(ga_models[i], 999)))
order_cd = sorted(range(N2), key=lambda i: (-_group_rank(ga_models[i]), -_alpha.get(ga_models[i], 999)))

fig2 = plt.figure(figsize=(MAX_WIDTH, MAX_HEIGHT), dpi=DPI, facecolor="white")
gs2 = gridspec.GridSpec(3, 2, height_ratios=[1, 0.14, 1],
    width_ratios=[1.15, 0.85], hspace=0.10, wspace=0.48,
    left=0.17, right=0.97, top=0.98, bottom=0.05)

ax_a  = fig2.add_subplot(gs2[0, 0]); ax_b  = fig2.add_subplot(gs2[0, 1])
ax_lg = fig2.add_subplot(gs2[1, :]); ax_lg.axis("off")
ax_c  = fig2.add_subplot(gs2[2, 0]); ax_d  = fig2.add_subplot(gs2[2, 1])

dumbbell(ax_a, ga_mcc,  amp_mcc, order_ab, "MCC", (-0.02, 0.88))
dumbbell(ax_b, ga_fpr,  amp_fpr, order_ab, "FPR", (-0.03, 1.05))
ax_a.set_title("AMP/non-AMP & GA", fontsize=FS_TITLE, pad=3)
ax_b.set_title("AMP/non-AMP & GA", fontsize=FS_TITLE, pad=3)
grouped_barh(ax_c, slay_p100, ga_p100, order_cd, "Precision@100", (0, 1.05), ref=slay_base)

sl = np.where(slay_lr > 0, slay_lr, 0.5)
gl = np.where(ga_lr   > 0, ga_lr,   0.5)
grouped_barh(ax_d, sl, gl, order_cd, "LR+ (log scale)", (0.4, 40), log=True, ref=1.0)
ax_c.set_title("GA & SLAY", fontsize=FS_TITLE, pad=3)
ax_d.set_title("GA & SLAY", fontsize=FS_TITLE, pad=3)
ax_d.set_xticks([0.5, 1, 2, 5, 10, 25])
ax_d.set_xticklabels(["0.5", "1", "2", "5", "10", "25"], fontsize=FS_TICK)
ax_d.xaxis.set_minor_formatter(mticker.NullFormatter())

_b = Patch(fc="none", ec="none", label=" ")
leg_shared = [
    # Group 1: model-type colors (4 items → fills one row)
    Patch(fc=CLF_CLR, ec="none", label="Classifier"),
    Patch(fc=REG_CLR, ec="none", label="Regressor"),
    Patch(fc=ACT_CLR, ec="none", label=r"HydrAMP$_{\mathrm{Ec}}$"),
    _b,
    # Group 2: bar fill style (4 items → fills one row)
    Patch(fc="#888", ec=GRID_CLR, lw=0.3, label="GA"),
    Patch(fc="#888", ec="white",  lw=0.3, hatch="//", label="SLAY"),
    _b, _b,
    # Group 3: dot styles — larger markers (4 items → fills one row)
    Line2D([0],[0], marker="o", color="w", mfc="#555", ms=7, label="GA"),
    Line2D([0],[0], marker="o", color="w", mfc="none", mec="#555", ms=7, mew=1.0,
           label="AMP / non-AMP"),
    _b, _b,
    # Group 4: median lines — three distinctive styles, no color (4 items → fills one row)
    Line2D([0],[0], color="#444", lw=1.2, ls="-",       alpha=0.7, label="Median (GA)"),
    Line2D([0],[0], color="#444", lw=1.2, ls=(0,(5,2)), alpha=0.7, label="Median (AMP/non-AMP)"),
    Line2D([0],[0], color="#444", lw=1.2, ls=(0,(1,1)), alpha=0.7, label="Median (SLAY)"),
    _b,
]
fig2.legend(handles=leg_shared, fontsize=FS_TICK, loc="center", ncol=4, frameon=False,
            bbox_to_anchor=(0.57, 0.504), bbox_transform=fig2.transFigure,
            handletextpad=0.3, columnspacing=0.9, handlelength=2.2, handleheight=0.7,
            borderpad=0.4, labelspacing=0.35)

for lbl, ax in [("a", ax_a), ("b", ax_b), ("c", ax_c), ("d", ax_d)]:
    p = ax.get_position()
    fig2.text(p.x0 - 0.07, p.y1 + 0.005, lbl,
              fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left")

save_figure(fig2, "figure2")
plt.show()
