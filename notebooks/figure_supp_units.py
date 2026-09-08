"""Supplementary figure: reporting MIC in ug/ml against uM.

MIC in ug/ml depends on molecular weight, so converting to uM reorders
peptides rather than rescaling them.  Panel a shows which peptides the unit
moves between activity classes, panel b what that does to model comparison,
and panels c and d what it does to the regression metrics.

All panels use the benchmark dataset (battleamp-all) and its two strain-level
MIC tasks, E. coli ATCC 25922 and S. aureus ATCC 25923.  The holdout-2026 and
deep-amp tasks are deliberately excluded.  Panel b's ug/ml labels reproduce the
published strain_ecoli25922 and strain_saureus25923 tasks to 99.9%, so it is
StrainActivity re-derived so that both units share one peptide set.

Panel b includes regressors binarised as classifiers, but they should be read
with care: changing the unit moves both their predictions and their labels,
so their shift reflects threshold placement as much as model quality.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

ROOT = Path(__file__).parent.parent / "results" / "unit_comparison"
reg = pd.read_csv(ROOT / "unit_metrics.tsv", sep="\t")
clf = pd.read_csv(ROOT / "classification_metrics.tsv", sep="\t")
pep = pd.read_csv(ROOT / "peptide_units.tsv", sep="\t")
eff = pd.read_csv(ROOT / "threshold_effect.tsv", sep="\t").set_index("task")

clf = clf[clf.variant != "example-model"]
SHOWN = "regression_ecoli25922"
# Two shades of one hue: the same measurement expressed two ways.  The
# active/inactive palette is reserved for activity categories.
UGML_COLOR, UM_COLOR = "#54278F", "#9E9AC8"
MARKERS = {"regression_ecoli25922": "o", "regression_saureus25923": "^"}
TASK_LABEL = {"regression_ecoli25922": "Ec ATCC 25922",
              "regression_saureus25923": "Sa ATCC 25923"}
TASK_SHORT = {"regression_ecoli25922": "Ec", "regression_saureus25923": "Sa"}

fig = plt.figure(figsize=(MAX_WIDTH, 6.1), dpi=DPI, facecolor="white")
gs = gridspec.GridSpec(2, 2, left=0.095, right=0.985, top=0.955, bottom=0.185,
                       hspace=0.34, wspace=0.26)


def unit_scatter(ax, table, metric, label, colour_by_type, title):
    piv = table.pivot_table(index=["variant", "task"], columns="unit",
                            values=metric, aggfunc="first")
    lo = hi = None
    for (variant, task), row in piv.iterrows():
        x, y = row.get("ug/ml"), row.get("uM")
        if pd.isna(x) or pd.isna(y):
            continue
        colour = model_color(variant) if colour_by_type else REG_CLR
        ax.plot(x, y, MARKERS[task], ms=4.6, color=colour, mec="white", mew=0.35,
                alpha=0.9, zorder=3)
        lo = x if lo is None else min(lo, x, y)
        hi = y if hi is None else max(hi, x, y)
    pad = 0.06 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", lw=0.6, ls="--", zorder=4)
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(f"{label}, " + r"$\mu$g/ml", fontsize=FS_LABEL)
    ax.set_ylabel(f"{label}, " + r"$\mu$M", fontsize=FS_LABEL)
    both = piv[["ug/ml", "uM"]].dropna()
    d = both["uM"] - both["ug/ml"]
    # Rank correlation per task: pooling the two strains mixes two orderings and
    # understates how stable the ranking is within either one.  r_s rather than
    # rho, because panel c already uses rho for the metric on its axes.
    ranks = []
    for task in MARKERS:
        rows = both[both.index.get_level_values("task") == task]
        if len(rows) > 2:
            ranks.append(f"{TASK_SHORT[task]} {rows['ug/ml'].corr(rows['uM'], method='spearman'):.2f}")
    ax.text(0.97, 0.04, f"median {d.median():+.3f}\n$r_s$  " + "  ".join(ranks),
            transform=ax.transAxes, fontsize=FS_TICK, va="bottom", ha="right",
            color="#333", linespacing=1.4)
    ax.set_title(title, fontsize=FS_TITLE, pad=4)
    style_ax(ax)


# a: which peptides change class, and how that depends on length
ax = fig.add_subplot(gs[0, 0])
sub = pep[pep.task == SHOWN].copy()
names = [r"$\leq$10", "11-20", "21-30", "31-40", ">40"]
sub["bin"] = pd.cut(sub.length, bins=[0, 10, 20, 30, 40, 200], labels=names, right=True)
active = {u: sub[f"mic_{u}"] <= 32 for u in ("ugml", "uM")}
only_ugml = (active["ugml"] & ~active["uM"]).groupby(sub["bin"], observed=False).mean() * 100
only_um = (~active["ugml"] & active["uM"]).groupby(sub["bin"], observed=False).mean() * 100
counts = sub.groupby("bin", observed=False).size()

x = np.arange(len(names))
ax.bar(x, only_ugml, width=0.7, color=UGML_COLOR, edgecolor="white", linewidth=0.3,
       zorder=3, label=r"active only in $\mu$g/ml")
ax.bar(x, only_um, width=0.7, bottom=only_ugml, color=UM_COLOR, edgecolor="white",
       linewidth=0.3, zorder=3, label=r"active only in $\mu$M")
for i, n in enumerate(counts):
    ax.text(i, only_ugml.iloc[i] + only_um.iloc[i] + 0.7, f"{n}", ha="center", va="bottom",
            fontsize=FS_TICK, color="#555")
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=FS_TICK)
ax.set_xlabel("peptide length (aa)", fontsize=FS_LABEL)
ax.set_ylabel("peptides changing\nactivity class (%)", fontsize=FS_LABEL)
ax.set_ylim(0, float(max(only_ugml + only_um)) * 1.3)
ax.set_title(f"Peptides reclassified, {TASK_LABEL[SHOWN]}", fontsize=FS_TITLE, pad=4)
ax.legend(fontsize=FS_TICK, frameon=False, loc="upper center", handlelength=1.2,
          handletextpad=0.4, borderpad=0.2)
style_ax(ax)

# b: what that does to model comparison
ax = fig.add_subplot(gs[0, 1])
unit_scatter(ax, clf, "mcc", "MCC", colour_by_type=True,
             title="StrainActivity")

# c, d: regression metrics
for j, (metric, label) in enumerate([("spearman", r"Spearman $\rho$"), ("r2_log2", r"$R^2_{\log_2}$")]):
    unit_scatter(fig.add_subplot(gs[1, j]), reg, metric, label, colour_by_type=False,
                 title="MIC regression")

# Two rows at the foot of the figure: strains, then model types.  Handle sizes
# follow figure2's legend strip.
strain_handles = [Line2D([0], [0], marker=MARKERS[t], color="w", mfc="#6e6e6e", mec="white",
                         ms=7, label=TASK_LABEL[t]) for t in MARKERS]
type_handles = [Patch(fc=c, ec="white", lw=0.3, label=l)
                for c, l in ((CLF_CLR, "Classifier"), (ACT_CLR, "Activity-aware"),
                             (REG_CLR, "Regressor"))]
for handles, ncol, y in ((strain_handles, 2, 0.083), (type_handles, 3, 0.045)):
    fig.legend(handles=handles, loc="center", ncol=ncol, fontsize=FS_TICK, frameon=False,
               bbox_to_anchor=(0.54, y), bbox_transform=fig.transFigure,
               handletextpad=0.3, columnspacing=0.9, handlelength=2.2, handleheight=0.7,
               borderpad=0.4, labelspacing=0.35)

for ax, letter in zip(fig.axes[:4], "abcd"):
    p = ax.get_position()
    fig.text(p.x0 - 0.062, p.y1 + 0.035, letter, fontsize=FS_PANEL, fontweight="bold", va="top")

save_figure(fig, "figure_supp_units")
