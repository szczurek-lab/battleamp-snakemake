"""Figure 4: What models have learned."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

# ── compute pipeline ──────────────────────────────────────────────────────────
RECOMPUTE_FIG4 = False  # Force recompute even if cache exists

SEQUENCE_FILES_FIG4 = {
    "AMP":             DATA_ROOT / "amp_positive.fasta",
    "GeneralActivity": DATA_ROOT / "activity" / "broad_positive.fasta",
    "Random":          DATA_ROOT / "syntax" / "synthetic_random.fasta",
    "Shuffled":        DATA_ROOT / "syntax" / "broad_positive_shuffled.fasta",
    "Realistic":       DATA_ROOT / "syntax" / "synthetic_realistic.fasta",
}

CDHIT_IDENTITY = 0.6; CDHIT_WORDSIZE = 4; MMSEQS_EVALUE = 1000
PHYSCHEM_SAMPLE_N = None; SIMILARITY_SAMPLE_N = 1000

LENGTH_TASKS  = ["length_01_10", "length_11_20", "length_21_30", "length_31_50"]
LENGTH_LABELS = ["1-10", "11-20", "21-30", "31-50"]
HOMOLOGY_TASKS = ["homology_40", "homology_60", "homology_80"]
HOM_LABELS     = ["40 %", "60 %", "80 %"]


def load_sequences(path, col="Sequence"):
    suffix = path.suffix.lower()
    if suffix in (".fasta", ".fa", ".faa"):
        seqs, seq = [], []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq: seqs.append("".join(seq)); seq = []
                else: seq.append(line)
            if seq: seqs.append("".join(seq))
        log.info(f"  Loaded {len(seqs):,} seqs from {path.name}")
        return seqs
    sep = "\t" if suffix == ".tsv" else ","
    df = pd.read_csv(path, sep=sep)
    return df[col].dropna().str.strip().tolist()


def compute_physchem(sequences, n=None):
    import modlamp.analysis as manalysis
    if n and n < len(sequences):
        sequences = list(np.random.choice(sequences, n, replace=False))
    props = {}
    for func, key in [
        (lambda s: (manalysis.GlobalAnalysis(s), "calc_charge", "charge"),  "charge"),
        (lambda s: (manalysis.GlobalDescriptor(s), "isoelectric_point", "descriptor"), "isoelectric_point"),
    ]:
        try:
            obj, method, attr = func(sequences)
            getattr(obj, method)()
            props[key] = list(getattr(obj, attr) if key == "charge" else obj.descriptor.flatten())
            if key == "charge": props[key] = list(obj.charge[0])
        except Exception as e:
            log.error(f"  {key}: {e}")
    try:
        h = manalysis.GlobalAnalysis(sequences); h.calc_H(scale="eisenberg")
        props["hydrophobicity"] = list(h.H[0])
    except Exception as e: log.error(f"  hydrophobicity: {e}")
    try:
        h = manalysis.PeptideDescriptor(sequences, "eisenberg"); h.calculate_moment()
        props["hydrophobic_moment"] = list(h.descriptor.flatten())
    except Exception as e: log.error(f"  hydrophobic_moment: {e}")
    try:
        h = manalysis.PeptideDescriptor(sequences, "levitt_alpha"); h.calculate_global()
        props["alpha_helix_propensity"] = list(h.descriptor.flatten())
    except Exception as e: log.error(f"  alpha_helix_propensity: {e}")
    return props


def load_predicted_positive_rates(model_ids, results_path=None, inference_dir=None, sequence_files=None):
    if results_path and results_path.exists():
        sep = "\t" if results_path.suffix == ".tsv" else ","
        df = pd.read_csv(results_path, sep=sep)
        if "variant" in df.columns: df = df.rename(columns={"variant": "model"})
        task_map = {"synthetic_random":"Random","synthetic_shuffled":"Shuffled","synthetic_realistic":"Realistic"}
        rates = {}
        for mid in model_ids:
            rates[mid] = {}
            for tf, ts in task_map.items():
                row = df[(df["model"]==mid) & (df["task"]==tf)]
                if row.empty: rates[mid][ts] = float("nan"); continue
                r = row.iloc[0]; nm = r.get("n_matched", r.get("n_samples", None)); np_ = r.get("n_predicted_positive", None)
                rates[mid][ts] = 100*np_/nm if pd.notna(nm) and pd.notna(np_) and nm>0 else float("nan")
        return rates
    log.warning("  Cannot load predicted positive rates"); return {}


def load_paired_scores(model_ids, positive_seqs, shuffled_seqs, inf_dir=INFERENCE_DIR):
    pos_set = set(positive_seqs); shuf_set = set(shuffled_seqs); results = {}
    for mid in model_ids:
        ip = inf_dir / mid / "battleamp-all" / "predictions.tsv"
        if not ip.exists(): results[mid] = {"positive":[],"shuffled":[]}; continue
        inf = _normalise_seq_col(pd.read_csv(ip, sep="\t"))
        if "sequence" not in inf.columns: results[mid]={"positive":[],"shuffled":[]}; continue
        if "Probability_score" in inf.columns: sc, is_c = "Probability_score", True
        elif "MIC" in inf.columns: sc, is_c = "MIC", False
        else: results[mid]={"positive":[],"shuffled":[]}; continue
        ps = pd.to_numeric(inf.loc[inf["sequence"].isin(pos_set), sc], errors="coerce").dropna().tolist()
        ss = pd.to_numeric(inf.loc[inf["sequence"].isin(shuf_set), sc], errors="coerce").dropna().tolist()
        if not is_c:
            ps = [np.log2(max(v,0.1)) for v in ps if v>0]
            ss = [np.log2(max(v,0.1)) for v in ss if v>0]
        results[mid] = {"positive": ps, "shuffled": ss}
    return results


def load_binned_mcc(model_ids, results_path=CLF_FILE):
    if not results_path.exists(): return {"length":{},"homology":{}}
    df = pd.read_csv(results_path, sep="\t")
    if "variant" in df.columns: df = df.rename(columns={"variant":"model"})
    lm, hm = {}, {}
    for mid in model_ids:
        lm[mid], hm[mid] = {}, {}
        for task in LENGTH_TASKS:
            row = df[(df["model"]==mid)&(df["task"]==task)]
            if row.empty: lm[mid][task]=float("nan"); continue
            r=row.iloc[0]; cov=r.get("coverage",1.0)
            lm[mid][task]=float("nan") if (pd.notna(cov) and float(cov)==0) or pd.isna(r.get("mcc")) else float(r["mcc"])
        for task in HOMOLOGY_TASKS:
            row = df[(df["model"]==mid)&(df["task"]==task)]
            if row.empty: hm[mid][task]=float("nan"); continue
            r=row.iloc[0]; cov=r.get("coverage",1.0)
            hm[mid][task]=float("nan") if (pd.notna(cov) and float(cov)==0) or pd.isna(r.get("mcc")) else float(r["mcc"])
    return {"length": lm, "homology": hm}


def compute_aa_frequencies(sequences, n=None):
    from collections import Counter
    if n and n<len(sequences): sequences=list(np.random.choice(sequences,n,replace=False))
    aa=Counter(); total=0
    for s in sequences: aa.update(s); total+=len(s)
    return {a:c/total for a,c in aa.items()} if total>0 else {}


def compute_and_cache_fig4():
    from metrics import compute_coverage, compute_similarity_to_reference
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    datasets = {}
    for name, path in SEQUENCE_FILES_FIG4.items():
        if path.exists(): datasets[name] = load_sequences(path)
        else: log.warning(f"  Not found: {path}")
    if not datasets: log.error("No datasets"); return

    log.info("Computing physicochemical properties...")
    physchem = {n: compute_physchem(s, n=PHYSCHEM_SAMPLE_N) for n, s in datasets.items()}
    json.dump(physchem, open(CACHE_DIR/"physchem.json","w"))

    log.info("Computing CD-HIT coverage...")
    cov = {n: compute_coverage(s, identity=CDHIT_IDENTITY, wordsize=CDHIT_WORDSIZE) for n,s in datasets.items()}
    json.dump(cov, open(CACHE_DIR/"coverage.json","w"))

    if "GeneralActivity" in datasets:
        log.info("Computing similarity to GeneralActivity...")
        ga = datasets["GeneralActivity"]
        sim = {n: compute_similarity_to_reference(s, ga, n=SIMILARITY_SAMPLE_N).tolist() for n,s in datasets.items()}
        json.dump(sim, open(CACHE_DIR/"sim_to_ga.json","w"))

    log.info("Computing AA frequencies...")
    json.dump({n: compute_aa_frequencies(s,1000) for n,s in datasets.items()}, open(CACHE_DIR/"aa_frequencies.json","w"))

    log.info("Loading predicted positive rates...")
    all_ids = list(CLF_FIG4.keys()) + list(REG_FIG4.keys())
    pp = load_predicted_positive_rates(all_ids, results_path=CLF_FILE if CLF_FILE.exists() else None,
                                       inference_dir=INFERENCE_DIR if INFERENCE_DIR.exists() else None,
                                       sequence_files=SEQUENCE_FILES_FIG4)
    json.dump(pp, open(CACHE_DIR/"predicted_positive.json","w"))

    log.info("Loading paired scores...")
    pos_s = datasets.get("GeneralActivity",[]); shuf_s = datasets.get("Shuffled",[])
    if pos_s and shuf_s:
        paired = load_paired_scores(list(MODEL_DISPLAY.keys()), pos_s, shuf_s)
        json.dump(paired, open(CACHE_DIR/"paired_scores.json","w"))

    log.info("Loading binned MCC...")
    binned = load_binned_mcc(HEATMAP_ORDER)
    json.dump(binned["length"], open(CACHE_DIR/"length_mcc.json","w"))
    json.dump(binned["homology"], open(CACHE_DIR/"homology_mcc.json","w"))
    log.info("Cache complete.")

_CACHE_FILES = ["physchem.json", "coverage.json", "sim_to_ga.json",
                "predicted_positive.json", "paired_scores.json",
                "length_mcc.json", "homology_mcc.json"]
_cache_missing = not all((CACHE_DIR / f).exists() for f in _CACHE_FILES)

if RECOMPUTE_FIG4 or _cache_missing:
    if _cache_missing:
        log.info("Figure 4 cache not found, computing from scratch...")
    else:
        log.info("Force-recomputing figure 4 cache...")
    compute_and_cache_fig4()
else:
    log.info("Figure 4 cache found, using existing files")

# ── plot from cache ───────────────────────────────────────────────────────────
def _lj(name):
    p = CACHE_DIR / name
    if p.exists(): return json.load(open(p))
    log.warning(f"  Missing: {p}"); return None

physchem = _lj("physchem.json") or {}
coverage = _lj("coverage.json") or {}
sim_to_ga = _lj("sim_to_ga.json") or {}
pp_rates = _lj("predicted_positive.json") or {}
paired_scores = _lj("paired_scores.json") or {}
length_mcc_d = _lj("length_mcc.json") or {}
homology_mcc_d = _lj("homology_mcc.json") or {}

ds_order = list(DATASET_PAL.keys())
ds_labels = ["AMP", "General\nActivity", "Random", "Shuffled", "Realistic"]
ds_colors = list(DATASET_PAL.values())
categories = ["Random", "Shuffled", "Realistic"]
BAR_COLORS_S = {"Random": RANDOM_COLOR, "Shuffled": SHUFFLED_COLOR, "Realistic": REALISTIC_COLOR}

def _sn4(mid): return CLF_FIG4.get(mid) or REG_FIG4.get(mid) or MODEL_DISPLAY.get(mid, mid)

def draw_violins(ax, data_list, colors, labels, ylabel, show_x=True):
    parts = ax.violinplot(data_list, showmeans=False, showmedians=False, showextrema=False,
                          widths=0.72, bw_method="silverman", points=80)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c); body.set_edgecolor("none"); body.set_alpha(0.85)
    for i, d in enumerate(data_list):
        q1,q2,q3 = np.percentile(d, [25,50,75]); x = i+1
        ax.vlines(x, q1, q3, color="black", lw=2.5, alpha=0.4, zorder=4)
        ax.scatter(x, q2, fc="white", ec="black", s=6, zorder=5, lw=0.4)
    ax.set_xticks(range(1, len(data_list)+1))
    ax.set_xticklabels(labels if show_x else [], rotation=55, ha="right", fontsize=FS_TICK)
    av = np.concatenate(data_list); lo,hi = np.quantile(av, 0.005), np.quantile(av, 0.995)
    ax.set_ylim(lo-(hi-lo)*0.12, hi+(hi-lo)*0.12)
    ax.set_ylabel(ylabel, fontsize=FS_LABEL); style_ax(ax)

fig4 = plt.figure(figsize=(MAX_WIDTH, MAX_HEIGHT), dpi=DPI)

# Vertical layout: b–c gap smaller (related panels); others standard
_Ha, _Hb, _Hc, _Hd = 0.092, 0.174, 0.092, 0.248
_gap_ab, _gap_bc, _gap_cd = 0.108, 0.061, 0.108
_a_top = 0.970;                _a_bot = _a_top - _Ha
_b_top = _a_bot - _gap_ab;    _b_bot = _b_top - _Hb
_c_top = _b_bot - _gap_bc;    _c_bot = _c_top - _Hc
_d_top = _c_bot - _gap_cd;    _d_bot = _d_top - _Hd
_LR = dict(left=0.09, right=0.88)

r1 = gridspec.GridSpec(1, 4, wspace=0.50, top=_a_top, bottom=_a_bot, **_LR)
for pi, (dn, key) in enumerate([("Charge","charge"),("Hydrophobic\nmoment","hydrophobic_moment"),
                                  ("Similarity\nto GA","_sim_to_ga")]):
    ax = fig4.add_subplot(r1[0, pi])
    data = []
    for ds in ds_order:
        vals = sim_to_ga.get(ds,[]) if key=="_sim_to_ga" else physchem.get(ds,{}).get(key,[])
        data.append(np.array(vals or [0.0], dtype=float))
    draw_violins(ax, data, ds_colors, ds_labels, dn)

axcov = fig4.add_subplot(r1[0, 3])
axcov.bar(range(len(ds_labels)), [coverage.get(ds,0) for ds in ds_order],
          color=ds_colors, ec="none", width=0.55, zorder=3)
axcov.set_xticks(range(len(ds_labels))); axcov.set_xticklabels(ds_labels, rotation=55, ha="right", fontsize=FS_TICK)
axcov.set_ylim(0, 1.08); axcov.yaxis.set_major_locator(plt.MultipleLocator(0.25))
axcov.set_ylabel("Diversity", fontsize=FS_LABEL); style_ax(axcov)

fig4.legend(handles=[Patch(fc=c, label=n) for n,c in DATASET_PAL.items()], loc="lower center",
    bbox_to_anchor=(0.48, _a_bot - 0.08), frameon=False,
    fontsize=FS_LEGEND, handlelength=1, handleheight=0.6,
    labelspacing=0.4, borderpad=0.3, ncol=5)
fig4.text(0.02, 0.97, "a", fontsize=FS_PANEL, fontweight="bold", va="top", transform=fig4.transFigure)

r2 = gridspec.GridSpec(1, 2, wspace=0.35, top=_b_top, bottom=_b_bot, width_ratios=[0.7, 1.3], **_LR)
def draw_hbars(ax, md, title):
    ns = list(md.values()); ids = list(md.keys()); nm = len(ns); y = np.arange(nm); h = 0.22
    for j, cat in enumerate(categories):
        vals = [pp_rates.get(mid,{}).get(cat,float("nan")) for mid in ids]
        ax.barh(y+(j-1)*h, vals, h, color=BAR_COLORS_S[cat], ec="none", zorder=3, label=cat)
    ax.set_yticks(y); ax.set_yticklabels(ns, fontsize=FS_TICK)
    ax.set_xlim(0, 105); ax.xaxis.set_major_locator(plt.MultipleLocator(25)); ax.invert_yaxis()
    ax.set_xlabel("Predicted positive (%)", fontsize=FS_LABEL); style_ax(ax)
    ax.grid(axis="x", color=GRID_CLR, lw=0.3, zorder=0); ax.grid(axis="y", visible=False)
    ax.set_title(title, fontsize=FS_TITLE, fontstyle="italic", pad=3)

draw_hbars(fig4.add_subplot(r2[0,0]), CLF_FIG4, "Classifiers")
draw_hbars(fig4.add_subplot(r2[0,1]), REG_FIG4, "Regressors")
fig4.text(0.02, _b_top + 0.012, "b", fontsize=FS_PANEL, fontweight="bold", va="top", transform=fig4.transFigure)

r3 = gridspec.GridSpec(1, 2, wspace=0.35, top=_c_top, bottom=_c_bot, width_ratios=[0.7, 1.3], **_LR)
def draw_paired(ax, md, title, is_clf=True):
    ids=list(md.keys()); ns=list(md.values()); data=[]; pos=[]; cols=[]; mt=[]
    for im, mid in enumerate(ids):
        ps=paired_scores.get(mid,{}); p=ps.get("positive",[]) or [0.0]; s=ps.get("shuffled",[]) or [0.0]
        b=im*3; data.extend([p,s]); pos.extend([b,b+1]); cols.extend([ACTIVE_COLOR, SHUFFLED_COLOR]); mt.append(b+0.5)
    bp=ax.boxplot(data, positions=pos, patch_artist=True, widths=0.75, showfliers=False,
        boxprops=dict(lw=0), medianprops=dict(lw=0.5,color="k"), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    for p,c in zip(bp["boxes"],cols): p.set_facecolor(c); p.set_alpha(0.9)
    ax.set_xticks(mt); ax.set_xticklabels(ns, rotation=40, ha="right", fontsize=FS_TICK)
    style_ax(ax); ax.grid(axis="x", visible=False)
    avl=[np.array(d,dtype=float) for d in data if len(d)>1]
    if avl:
        av=np.concatenate(avl); lo,hi=np.quantile(av,0.005),np.quantile(av,0.995); pad=(hi-lo)*0.08
        ax.set_ylim(lo-pad, hi+pad)
    ax.set_ylabel("$P_{AMP}$" if is_clf else r"$\log_2$(MIC)", fontsize=FS_LABEL)
    ax.set_title(title, fontsize=FS_TITLE, fontstyle="italic", pad=3)

ax_cp = fig4.add_subplot(r3[0,0])
draw_paired(ax_cp, CLF_FIG4, "Classifiers", True)
ax_cp.legend(handles=[Patch(fc=ACTIVE_COLOR,label="Positive"),Patch(fc=SHUFFLED_COLOR,label="Shuffled")],
    loc="lower left", frameon=True, fancybox=False, edgecolor="0.7", framealpha=0.9,
    fontsize=4.5, handlelength=1, handleheight=0.6)
ax_rp = fig4.add_subplot(r3[0,1])
draw_paired(ax_rp, REG_FIG4, "Regressors", False)
ax_rp.axhline(y=np.log2(32), color="#aaa", lw=0.9, ls="--", alpha=0.9, zorder=0)
ax_rp.legend(handles=[Line2D([0],[0], color="#aaa", lw=0.9, ls="--", label=r"MIC = 32 $\mu$g/ml")],
    loc="upper left", frameon=True, fancybox=False, edgecolor="0.7", framealpha=0.9,
    fontsize=FS_TICK, handlelength=1.5, handleheight=0.6, borderpad=0.3)
fig4.text(0.02, _c_top + 0.012, "c", fontsize=FS_PANEL, fontweight="bold", va="top", transform=fig4.transFigure)

r4 = gridspec.GridSpec(1, 2, wspace=0.12, top=_d_top, bottom=_d_bot, width_ratios=[1, 0.8], left=0.18, right=0.88)
nhm = len(HEATMAP_ORDER); ncl = len(CLF_FIG4)
lm = np.full((nhm, len(LENGTH_TASKS)), np.nan); hm = np.full((nhm, len(HOMOLOGY_TASKS)), np.nan)
for i, mid in enumerate(HEATMAP_ORDER):
    for j, t in enumerate(LENGTH_TASKS):
        if mid in length_mcc_d: lm[i,j] = length_mcc_d[mid].get(t, np.nan)
    for j, t in enumerate(HOMOLOGY_TASKS):
        if mid in homology_mcc_d: hm[i,j] = homology_mcc_d[mid].get(t, np.nan)

median_lm_row = np.nanmedian(lm, axis=0, keepdims=True)
median_hm_row = np.nanmedian(hm, axis=0, keepdims=True)
# gap1: classifiers/regressors; gap2: before aggregated rows
_gap_l = np.full((1, len(LENGTH_TASKS)),   np.nan)
_gap_h = np.full((1, len(HOMOLOGY_TASKS)), np.nan)
lm_ext = np.vstack([lm[:ncl], _gap_l, lm[ncl:], _gap_l.copy(), median_lm_row])
hm_ext = np.vstack([hm[:ncl], _gap_h, hm[ncl:], _gap_h.copy(), median_hm_row])
nhm_ext  = lm_ext.shape[0]   # nhm + 3
GAP2     = nhm + 1
_GAP_ROWS = {ncl, GAP2}
_AGG_ROWS = {nhm + 2}
yticklabels_d = ([_sn4(m) for m in HEATMAP_ORDER[:ncl]] + [""] +
                 [_sn4(m) for m in HEATMAP_ORDER[ncl:]] + ["", "Median"])
_model_rows = [i for i in range(nhm_ext) if i not in _GAP_ROWS and i not in _AGG_ROWS]

cm4 = matplotlib.colormaps["RdYlGn"]; cm4.set_bad("#e8e8e8")
nm4 = mcolors.TwoSlopeNorm(vmin=-0.1, vcenter=0.3, vmax=0.7)

# Figure-coordinate y-centres for group labels (row i centred at data-y=i, y-axis inverted)
_row_h_d    = (_d_top - _d_bot) / nhm_ext
_clf_ctr_fig = _d_top - (ncl / 2) * _row_h_d
_reg_ctr_fig = _d_top - ((ncl + nhm + 2) / 2) * _row_h_d

for mi, (mat, labs, title, show_y) in enumerate([
    (lm_ext, LENGTH_LABELS, "MCC by peptide length", True),
    (hm_ext, HOM_LABELS, "MCC by max homology to training data", False)]):
    ax = fig4.add_subplot(r4[0, mi])
    for sp in ax.spines.values(): sp.set_visible(False)
    im = ax.imshow(mat, aspect="auto", cmap=cm4, norm=nm4, interpolation="nearest")
    for _gr in _GAP_ROWS:
        ax.add_patch(mpatches.Rectangle((-0.5, _gr - 0.5), mat.shape[1], 1.0,
                     fc="white", ec="none", zorder=3, transform=ax.transData, clip_on=True))
    for i in range(mat.shape[0]):
        if i in _GAP_ROWS: continue
        for j in range(mat.shape[1]):
            v = mat[i,j]
            if np.isnan(v): ax.text(j,i,"-",ha="center",va="center",fontsize=FS_TICK,color="#999",zorder=4)
            else: ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=FS_TICK,
                         color="white" if v<0.1 else "black", zorder=4)
    for j in range(mat.shape[1]):
        col = [(mat[i, j], i) for i in _model_rows if not np.isnan(mat[i, j])]
        if col:
            best = max(v for v, _ in col)
            for v, bi in col:
                if np.isclose(v, best):
                    ax.add_patch(mpatches.Rectangle((j - 0.5, bi - 0.5), 1, 1,
                                 fc="none", ec="black", lw=1.0, zorder=5))
    ax.set_xticks(np.arange(len(labs))); ax.set_xticklabels(labs, fontsize=FS_TICK)
    ax.set_yticks(np.arange(nhm_ext))
    ax.set_yticklabels(yticklabels_d if show_y else [], fontsize=FS_TICK)
    if show_y:
        tl = ax.get_yticklabels()
        tl[-1].set_fontweight("bold"); tl[-1].set_color("#555")
    ax.tick_params(axis="both", length=0); ax.set_title(title, fontsize=FS_TITLE, pad=4, color="#555")

fig4.text(0.02, _clf_ctr_fig, "Classifiers", ha="center", va="center",
          fontsize=FS_TICK, rotation=90, transform=fig4.transFigure, clip_on=False)
fig4.text(0.02, _reg_ctr_fig, "Regressors\n(binarized)", ha="center", va="center",
          fontsize=FS_TICK, rotation=90, transform=fig4.transFigure, clip_on=False)

cba = fig4.add_axes([0.90, _d_bot, 0.012, _d_top - _d_bot])
cb4 = fig4.colorbar(im, cax=cba); cb4.ax.tick_params(labelsize=FS_TICK); cb4.set_label("MCC", fontsize=FS_LABEL)
fig4.text(0.02, _d_top + 0.012, "d", fontsize=FS_PANEL, fontweight="bold", va="top", transform=fig4.transFigure)

save_figure(fig4, "figure4")
plt.show()
