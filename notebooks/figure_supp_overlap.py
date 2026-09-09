"""Supplementary figure: training-data overlap and its effect on GeneralActivity.

Three evaluation sets, one task (broad_activity / GeneralActivity):

  BattleAMP          battleamp-all, 4,355 sequences
  BattleAMP - train  the same task with each model's own training sequences
                     removed, so the set size differs per model (1,977-4,186)
  Holdout-2026       DBAASP peptides curated after every model was published,
                     521 sequences, seen by none of them

Overlap with each model's training data is given as a percentage beside its name.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

# Built from results/generalization/generalization_metrics.tsv, which is tracked,
# so this figure regenerates from a clone. Pass a pre-pivoted table as argv[1] to
# override. The source is long format (one row per variant, task and evaluation
# set); this pivots the broad_activity rows into the wide layout the panels want:
# one row per variant, columns <set>_<metric>, plus overlap_pct.
GEN_FILE = RESULTS_ROOT / "generalization" / "generalization_metrics.tsv"


def _build_overlap_table():
    g = pd.read_csv(GEN_FILE, sep="\t")
    g = g[g["task"] == "broad_activity"].copy()
    # LR+ = TPR/FPR. FPR of zero means no false positives, an infinite ratio;
    # the panel caps it at LR_CAP and flags the bar, so inf is the honest value.
    g["lr"] = np.where(g["fpr"] > 0, g["tpr"] / g["fpr"], np.inf)
    wide = g.pivot(index="variant", columns="evaluation_set",
                   values=["mcc", "lr", "fpr"])
    wide.columns = [f"{s}_{metric}" for metric, s in wide.columns]
    # overlap_pct is a property of the variant, identical across sets.
    wide["overlap_pct"] = g.groupby("variant")["overlap_pct"].first()
    return wide


if len(sys.argv) > 1:
    df = pd.read_csv(Path(sys.argv[1]), sep="\t").set_index("variant")
else:
    df = _build_overlap_table()
order = [m for m in ALL_MODELS if m in df.index]

SETS = [("battleamp", "BattleAMP", 1.00),
        ("minustrain", "BattleAMP $-$ train", 0.62),
        ("holdout", "Holdout-2026", 0.32)]
PANELS = [("mcc", "MCC", "linear"), ("lr", "LR$+$", "log"), ("fpr", "FPR", "linear")]
LR_CAP = 40.0          # display cap; LR+ is unbounded when FPR reaches 0

fig = plt.figure(figsize=(MAX_WIDTH, 5.6), dpi=DPI, facecolor="white")
gs = gridspec.GridSpec(1, 3, left=0.235, right=0.99, top=0.955, bottom=0.135,
                       wspace=0.14)

h = 0.27
ypos = np.arange(len(order))
axes = []

axa = None
for j, (metric, label, scale) in enumerate(PANELS):
    ax = fig.add_subplot(gs[0, j], sharey=axa) if axa is not None else fig.add_subplot(gs[0, j])
    if axa is None:
        axa = ax
    axes.append(ax)
    for k, (prefix, _, alpha) in enumerate(SETS):
        vals, capped = [], []
        for m in order:
            v = df.loc[m, f"{prefix}_{metric}"]
            capped.append(metric == "lr" and (v == np.inf or v > LR_CAP))
            vals.append(min(v, LR_CAP) if metric == "lr" else v)
        base = 1.0 if scale == "log" else 0.0
        ax.barh(ypos + (k - 1) * h, np.array(vals) - base, height=h, left=base,
                color=[model_color(m) for m in order], alpha=alpha,
                edgecolor="white", linewidth=0.3, zorder=3)
        for i, c in enumerate(capped):
            if c:
                ax.text(LR_CAP, ypos[i] + (k - 1) * h, r"$\infty$", ha="left",
                        va="center", fontsize=FS_TICK, color="#444", zorder=5)
    ax.set_xlabel(label, fontsize=FS_LABEL)
    if scale == "log":
        ax.set_xscale("log"); ax.set_xlim(0.8, LR_CAP * 1.5)
        ax.set_xticks([1, 3, 10, 30]); ax.set_xticklabels(["1", "3", "10", "30"])
        ax.axvline(1, color="black", lw=0.6, zorder=4)   # LR+ = 1 is uninformative
    elif metric == "fpr":
        ax.set_xlim(0, 1.02); ax.set_xticks([0, 0.5, 1.0])
    else:
        ax.axvline(0, color="black", lw=0.6, zorder=4)
    if ax is not axa:
        plt.setp(ax.get_yticklabels(), visible=False)

for ax in axes:
    ax.set_facecolor(BG)
    ax.tick_params(axis="both", length=0, labelsize=FS_TICK)
    ax.grid(axis="x", color=GRID_CLR, lw=0.3)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    for i in range(0, len(order), 2):
        ax.axhspan(i - 0.5, i + 0.5, color="white", alpha=0.55, zorder=0)

axa.set_yticks(ypos)
axa.set_yticklabels([f"{short_name(m)}  ({df.loc[m, 'overlap_pct']:.0f}%)" for m in order],
                    fontsize=FS_TICK)
for lbl, m in zip(axa.get_yticklabels(), order):
    lbl.set_color(model_color(m))
axa.set_ylim(len(order) - 0.5, -0.5)

for ax, letter in zip(axes, "abc"):
    p = ax.get_position()
    fig.text(p.x0 - 0.035, p.y1 + 0.012, letter, fontsize=FS_PANEL, fontweight="bold", va="bottom")

shade = [Patch(facecolor="#6e6e6e", alpha=a, edgecolor="white", linewidth=0.3, label=l)
         for _, l, a in SETS]
kind = [Patch(facecolor=c, edgecolor="white", linewidth=0.3, label=l)
        for c, l in ((CLF_CLR, "Classifier"), (ACT_CLR, "Activity-aware"), (REG_CLR, "Regressor"))]
# Centred on the subplot block (left..right), not on the figure, so it sits
# under the panels rather than under the panels plus the label gutter.
_p0, _p1 = axes[0].get_position(), axes[-1].get_position()
fig.legend(handles=shade + kind, loc="upper center", ncol=6, fontsize=FS_TICK,
           frameon=False, handlelength=1.3, handletextpad=0.4, columnspacing=1.4,
           bbox_to_anchor=((_p0.x0 + _p1.x1) / 2, 0.052))

save_figure(fig, "figure_supp_overlap")
