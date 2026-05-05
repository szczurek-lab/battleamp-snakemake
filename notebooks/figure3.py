"""Figure 3: Biological hierarchy evaluation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False

CLF_TASKS = [
    "broad_activity", "gram_minus", "gram_plus",
    "species_ecoli", "species_saureus", "species_paeruginosa",
    "species_kpneumoniae", "species_abaumannii",
    "strain_ecoli25922", "strain_saureus25923", "strain_saureus33591",
    "strain_saureus43300", "strain_paeruginosa27853",
    "strain_kpneumoniae700603", "strain_abaumannii19606",
]
CLF_TASK_SHORT = [
    "General\nActivity", "Gram$-$", "Gram$+$",
    "Ec", "Sa", "Pa", "Kp", "Ab",
    "Ec\n25922", "Sa\n25923", "Sa\n33591", "Sa\n43300",
    "Pa\n27853", "Kp\n700603", "Ab\n19606",
]
TASK_GROUPS = [("General", 0, 1), ("Gram", 1, 3), ("Species", 3, 8), ("Strain", 8, 15)]
N_CLF = len(CLASSIFIERS); N_ACT = len(ACTIVITY_AWARE)
MODEL_GROUPS = [("Classifiers", 0, N_CLF + N_ACT),
                ("Regressors (binarized)", N_CLF + N_ACT, len(ALL_MODELS))]
TRAINING_TARGET = {
    "sensexamp-ecoli": [3], "sensexamp-saureus": [4],
    "deep-amp-lstm-gramneg": [1], "deep-amp-lstm-grampos": [2],
    "deep-amp-cnn-gramneg": [1], "deep-amp-cnn-grampos": [2],
    "apex-ecoli": [3], "apex-saureus": [4], "apex-abaumannii": [7],
    "apex-paeruginosa": [5], "apex-kpneumoniae": [6],
    "mbc-attention": [3, 4], "ampredictor": [3], "hydramp-mic-classifier": [3, 4],
}

REG_TASKS = ["regression_ecoli25922", "regression_saureus25923"]
REG_TASK_LABELS = ["Ec ATCC 25922", "Sa ATCC 25923"]
REG_METRICS = ["r2_log2", "spearman", "msl2e"]
REG_METRIC_LABELS = [r"$R^2_{\log_2}$", r"$\rho$", "MSL2E"]
HIGHER_BETTER = [True, True, False]

df_clf = pd.read_csv(CLF_FILE, sep="\t"); df_clf = df_clf[df_clf["variant"] != "example-model"].copy()
df_reg = pd.read_csv(REG_FILE, sep="\t"); df_reg = df_reg[df_reg["variant"] != "example-model"].copy()
df_all = pd.concat([df_clf, df_reg], ignore_index=True)

def build_matrix(models, tasks, metric, data=None):
    if data is None: data = df_all
    mat = np.full((len(models), len(tasks)), np.nan)
    for i, m in enumerate(models):
        for j, t in enumerate(tasks):
            row = data[(data["variant"] == m) & (data["task"] == t)]
            if not row.empty:
                val = row.iloc[0][metric]
                if pd.notna(val) and val != "":
                    try: mat[i, j] = float(val)
                    except (ValueError, TypeError): pass
    return mat

clf_mcc = build_matrix(ALL_MODELS, CLF_TASKS, "mcc")
clf_cov = build_matrix(ALL_MODELS, CLF_TASKS, "coverage")

NCL = N_CLF + N_ACT   # 7: boundary between classifiers and regressors

fig3 = plt.figure(figsize=(MAX_WIDTH, MAX_HEIGHT), dpi=DPI, facecolor="white")
n_tasks = len(CLF_TASKS); n_models = len(ALL_MODELS)
N_CLF_ROWS = NCL; N_REG_ROWS = n_models - NCL

# Two subplots for panel a: classifier block / half-row gap / regressor block
_row_h  = (0.97 - 0.56) / (N_CLF_ROWS + 0.5 + N_REG_ROWS)
_clf_bot = 0.97 - N_CLF_ROWS * _row_h
_reg_top = _clf_bot - 0.5 * _row_h
gs_clf = gridspec.GridSpec(1, 1, left=0.14, right=0.92, top=0.97,      bottom=_clf_bot)
gs_reg = gridspec.GridSpec(1, 1, left=0.14, right=0.92, top=_reg_top,  bottom=0.56)
gs_b3  = gridspec.GridSpec(1, 1, left=0.19, right=0.87, top=0.48,      bottom=0.31)
gs_c3  = gridspec.GridSpec(1, 2, left=0.10, right=0.95, top=0.24,      bottom=0.03, wspace=0.35)

ax3a_clf = fig3.add_subplot(gs_clf[0, 0])
ax3a_reg = fig3.add_subplot(gs_reg[0, 0])

mcc_clf = clf_mcc[:NCL, :]; cov_clf = clf_cov[:NCL, :]
mcc_reg = clf_mcc[NCL:, :]; cov_reg = clf_cov[NCL:, :]
mask_clf = np.isnan(mcc_clf) | ((~np.isnan(cov_clf)) & (cov_clf < 0.5))
mask_reg = np.isnan(mcc_reg) | ((~np.isnan(cov_reg)) & (cov_reg < 0.5))
norm_a = TwoSlopeNorm(vmin=-0.2, vcenter=0, vmax=0.75)

sns.heatmap(mcc_clf, mask=mask_clf, ax=ax3a_clf, cmap="RdYlGn", norm=norm_a,
    linewidths=0.2, linecolor="white", annot=True, fmt=".2f",
    annot_kws={"size": FS_ANNOT}, cbar=False,
    xticklabels=False, yticklabels=[short_name(m) for m in ALL_MODELS[:NCL]])
sns.heatmap(mcc_reg, mask=mask_reg, ax=ax3a_reg, cmap="RdYlGn", norm=norm_a,
    linewidths=0.2, linecolor="white", annot=True, fmt=".2f",
    annot_kws={"size": FS_ANNOT}, cbar=False,
    xticklabels=CLF_TASK_SHORT, yticklabels=[short_name(m) for m in ALL_MODELS[NCL:]])

for i in range(N_CLF_ROWS):
    for j in range(n_tasks):
        if mask_clf[i, j]:
            ax3a_clf.add_patch(mpatches.Rectangle((j, i), 1, 1, fc="#f0f0f0", ec="white", lw=0.2))
for i in range(N_REG_ROWS):
    for j in range(n_tasks):
        if mask_reg[i, j]:
            ax3a_reg.add_patch(mpatches.Rectangle((j, i), 1, 1, fc="#f0f0f0", ec="white", lw=0.2))

for mk, ti in TRAINING_TARGET.items():
    if mk in ALL_MODELS:
        mi = ALL_MODELS.index(mk)
        if mi < NCL:
            ax_t, local_i, mask_t = ax3a_clf, mi, mask_clf
        else:
            ax_t, local_i, mask_t = ax3a_reg, mi - NCL, mask_reg
        for j in ti:
            if not mask_t[local_i, j]:
                ax_t.plot(j + 0.85, local_i + 0.15, "o", ms=2.5, color="black", mec="white", mew=0.3, zorder=5)

for j in range(n_tasks):
    cv_c = mcc_clf[:, j].copy(); cv_c[mask_clf[:, j]] = np.nan
    cv_r = mcc_reg[:, j].copy(); cv_r[mask_reg[:, j]] = np.nan
    bc = np.nanmax(cv_c) if not np.all(np.isnan(cv_c)) else -np.inf
    br = np.nanmax(cv_r) if not np.all(np.isnan(cv_r)) else -np.inf
    if bc == br == -np.inf: continue
    if bc >= br:
        ax3a_clf.add_patch(mpatches.Rectangle((j, np.nanargmax(cv_c)), 1, 1, fc="none", ec="black", lw=1.2, zorder=4))
    else:
        ax3a_reg.add_patch(mpatches.Rectangle((j, np.nanargmax(cv_r)), 1, 1, fc="none", ec="black", lw=1.2, zorder=4))

# Gram stain colours — intermediate between dark and fully saturated
_GRAM_NEG = "#1565C0"   # medium-dark blue (Gram-)
_GRAM_POS = "#C62828"   # medium-dark red  (Gram+)
# index → colour (None = neutral); matches CLF_TASKS order
_GRAM_CLR = [None, _GRAM_NEG, _GRAM_POS,
             _GRAM_NEG, _GRAM_POS, _GRAM_NEG, _GRAM_NEG, _GRAM_NEG,
             _GRAM_NEG, _GRAM_POS, _GRAM_POS, _GRAM_POS,
             _GRAM_NEG, _GRAM_NEG, _GRAM_NEG]
_MDR_IDX = {10, 11}  # Sa ATCC 33591 and 43300 — confirmed MRSA in dataset

ax3a_reg.set_xticklabels(CLF_TASK_SHORT, rotation=45, ha="right", fontsize=FS_TICK)
ax3a_reg.xaxis.set_ticks_position("bottom")
for j, lbl in enumerate(ax3a_reg.get_xticklabels()):
    if _GRAM_CLR[j] is not None:
        lbl.set_color(_GRAM_CLR[j])
    if j in _MDR_IDX:
        lbl.set_fontweight("bold")
for ax in [ax3a_clf, ax3a_reg]:
    for lbl in ax.get_yticklabels(): lbl.set_fontsize(FS_TICK); lbl.set_color("black")
    ax.tick_params(axis="both", length=0)

_cax  = fig3.add_axes([0.930, 0.56, 0.012, 0.97 - 0.56])
_sm   = plt.cm.ScalarMappable(cmap="RdYlGn", norm=norm_a); _sm.set_array([])
_cbar = fig3.colorbar(_sm, cax=_cax)
_cbar.ax.tick_params(labelsize=FS_TICK); _cbar.set_label("MCC", fontsize=FS_LABEL)

for ax in [ax3a_clf, ax3a_reg]:
    for _, _, e in TASK_GROUPS[:-1]: ax.axvline(x=e, color="black", lw=0.5, alpha=0.6)
for lb, s, e in TASK_GROUPS:
    ax3a_clf.text((s+e)/2/n_tasks, 1.01, lb, ha="center", va="bottom", fontsize=FS_TICK, transform=ax3a_clf.transAxes)
ax3a_clf.text(-0.18, 0.5, "Classifiers", ha="center", va="center",
              fontsize=FS_TICK, rotation=90, transform=ax3a_clf.transAxes, clip_on=False)
ax3a_reg.text(-0.18, 0.5, "Regressors (binarized)", ha="center", va="center",
              fontsize=FS_TICK, rotation=90, transform=ax3a_reg.transAxes, clip_on=False)

ax3b = fig3.add_subplot(gs_b3[0, 0])
n_orgs = len(REG_TASKS); n_mets = len(REG_METRICS)
reg_raw = np.full((len(REGRESSORS), n_orgs * n_mets), np.nan)
for i, m in enumerate(REGRESSORS):
    for j, t in enumerate(REG_TASKS):
        for k, met in enumerate(REG_METRICS):
            row = df_reg[(df_reg["variant"] == m) & (df_reg["task"] == t)]
            if not row.empty:
                val = row.iloc[0][met]
                if pd.notna(val) and val != "":
                    try: reg_raw[i, j*n_mets+k] = float(val)
                    except: pass

reg_norm = reg_raw.copy()
for ci in range(reg_raw.shape[1]):
    mi2 = ci % n_mets; col = reg_raw[:, ci]; v = col[~np.isnan(col)]
    if len(v) > 0 and v.max() > v.min():
        reg_norm[:, ci] = (col - v.min()) / (v.max() - v.min())
        if not HIGHER_BETTER[mi2]: reg_norm[:, ci] = 1 - reg_norm[:, ci]
    elif len(v) > 0: reg_norm[:, ci] = 0.5

at = np.full(reg_raw.shape, "", dtype=object)
for i in range(reg_raw.shape[0]):
    for j in range(reg_raw.shape[1]):
        if not np.isnan(reg_raw[i,j]):
            at[i,j] = f"{reg_raw[i,j]:.1f}" if j%n_mets==2 else f"{reg_raw[i,j]:.2f}"

mask_b = np.isnan(reg_raw)
cl = [ml for _ in REG_TASK_LABELS for ml in REG_METRIC_LABELS]
sns.heatmap(reg_norm, mask=mask_b, ax=ax3b, cmap="RdYlGn", vmin=0, vmax=1,
    linewidths=0.2, linecolor="white", annot=at, fmt="", annot_kws={"size": FS_ANNOT},
    cbar=False, xticklabels=cl, yticklabels=[short_name(m) for m in REGRESSORS])
ax3b.set_xticklabels(ax3b.get_xticklabels(), rotation=0, ha="center", fontsize=FS_TICK)
for lbl in ax3b.get_yticklabels(): lbl.set_fontsize(FS_TICK); lbl.set_color("black")
ax3b.tick_params(axis="both", length=0)
for i in range(len(REGRESSORS)):
    for j in range(n_orgs*n_mets):
        if mask_b[i,j]: ax3b.add_patch(mpatches.Rectangle((j,i),1,1,fc="#f0f0f0",ec="white",lw=0.2))
for j in range(n_orgs*n_mets):
    col = [(reg_norm[i, j], at[i, j], i) for i in range(len(REGRESSORS)) if not mask_b[i, j] and at[i, j] != '']
    if col:
        _, best_ann, _ = max(col, key=lambda x: x[0])
        for _, ann, bi in col:
            if ann == best_ann:
                ax3b.add_patch(mpatches.Rectangle((j, bi), 1, 1, fc="none", ec="black", lw=1.0, zorder=5))
ax3b.axvline(x=n_mets, color="black", lw=1)
for idx, ol in enumerate(REG_TASK_LABELS):
    _tc = _GRAM_NEG if idx == 0 else _GRAM_POS  # Ec=Gram-, Sa=Gram+
    ax3b.text((idx*n_mets+n_mets/2)/(n_mets*n_orgs), 1.04, ol, ha="center", va="bottom",
              fontsize=FS_LABEL, transform=ax3b.transAxes, color=_tc)
bp = ax3b.get_position()
cax_b = fig3.add_axes([bp.x1+0.02, bp.y0, 0.012, bp.height])
sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=plt.Normalize(0,1)); sm.set_array([])
cb = fig3.colorbar(sm, cax=cax_b); cb.set_ticks([0,0.5,1])
cax_b.set_yticklabels(["Worst","","Best"]); cax_b.tick_params(labelsize=FS_TICK, size=2)
cb.set_label("Per-column rank", fontsize=FS_TICK, labelpad=2)

rd = []
for m in REGRESSORS:
    for ti, t in enumerate(REG_TASKS):
        row = df_reg[(df_reg["variant"]==m)&(df_reg["task"]==t)]
        if not row.empty:
            r2 = row.iloc[0].get("r2_log2", np.nan); rho = row.iloc[0].get("spearman", np.nan)
            if pd.notna(r2) and pd.notna(rho):
                rd.append({"model":m,"task":t,"task_idx":ti,"r2":float(r2),"rho":float(rho)})
rsd = pd.DataFrame(rd)
EC_TR = {"mbc-attention","ampredictor","sensexamp-ecoli","apex-ecoli"}
SA_TR = {"mbc-attention","sensexamp-saureus","apex-saureus"}

axc = []
for ti, (tk, tl) in enumerate(zip(REG_TASKS, REG_TASK_LABELS)):
    ax = fig3.add_subplot(gs_c3[0, ti]); axc.append(ax)
    ax.set_facecolor(BG)
    sub = rsd[rsd["task"]==tk]; ts = EC_TR if ti==0 else SA_TR
    texts=[]; xc=[]; yc=[]
    for _, r in sub.iterrows():
        tr = r["model"] in ts; mk = "D" if tr else "o"; ms = 4.5 if tr else 3.5
        c = "#2166ac" if tr else "#b2182b"
        ax.plot(r["rho"], r["r2"], marker=mk, ms=ms, color=c, mec="white", mew=0.4, alpha=0.9, zorder=3)
        xc.append(r["rho"]); yc.append(r["r2"])
        texts.append(ax.text(r["rho"], r["r2"], scatter_name(r["model"]), fontsize=FS_TICK, ha="left", va="center", zorder=4))
    ax.axhline(y=0, color="black", lw=0.6, ls="--", alpha=0.5)
    rr = sub["rho"].values; r2r = sub["r2"].values
    ax.set_xlim(max(-0.08, rr.min()-max((rr.max()-rr.min())*0.6, 0.12)),
                min(1.10, rr.max()+max((rr.max()-rr.min())*0.6, 0.12)))
    ax.set_ylim(r2r.min()-max((r2r.max()-r2r.min())*0.4, 0.20),
                r2r.max()+max((r2r.max()-r2r.min())*0.4, 0.20))
    if HAS_ADJUST:
        adjust_text(texts, x=xc, y=yc, ax=ax, arrowprops=dict(arrowstyle="-", color="#999", lw=0.4),
            force_text=(1.2,1.2), force_static=(0.8,0.8), force_pull=(0.01,0.01),
            force_explode=(0.6,0.6), expand=(2.0,2.5), max_move=(10,10),
            iter_lim=5000, time_lim=30,
            ensure_inside_axes=True, prevent_crossings=True)
    ax.set_xlabel(r"Spearman $\rho$", fontsize=FS_LABEL)
    if ti==0: ax.set_ylabel(r"$R^2_{\log_2}$", fontsize=FS_LABEL)
    _tc = _GRAM_NEG if ti == 0 else _GRAM_POS  # Ec=Gram-, Sa=Gram+
    ax.set_title(tl, fontsize=FS_LABEL, pad=4, color=_tc)
    style_ax(ax)

axc[1].legend(handles=[
    Line2D([0],[0],marker="D",color="w",mfc="#2166ac",mec="white",ms=5,label="Target-matched"),
    Line2D([0],[0],marker="o",color="w",mfc="#b2182b",mec="white",ms=4,label="Cross-target"),
    Line2D([0],[0],color="black",lw=0.6,ls="--",label="$R^2 = 0$"),
], loc="lower right", fontsize=FS_TICK, frameon=True, facecolor="white", edgecolor="#ccc",
   handlelength=1.5, handletextpad=0.4)

fig3.text(0.02, 0.975, "a", fontsize=FS_PANEL, fontweight="bold", va="top")
fig3.text(0.05, 0.495, "b", fontsize=FS_PANEL, fontweight="bold", va="top")
fig3.text(0.02, 0.255, "c", fontsize=FS_PANEL, fontweight="bold", va="top")
save_figure(fig3, "figure3")
plt.show()
