"""Supplementary: True vs predicted MIC scatter plots."""
import sys
import math
from pathlib import Path
from scipy.stats import spearmanr
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

MIC_CACHE = PROJECT_ROOT / "figure_supp_mic_cache"

STRAIN_TRUTH = OrderedDict([
    ("ecoli", {
        "path": DATA_ROOT / "activity" / "strain" / "escherichiacoliatcc25922_mic.csv",
        "label": "E. coli ATCC 25922", "suffix": "ecoli",
    }),
    ("saureus", {
        "path": DATA_ROOT / "activity" / "strain" / "staphylococcusaureusatcc25923_mic.csv",
        "label": "S. aureus ATCC 25923", "suffix": "saureus",
    }),
])

def _find_mic_col(df, model_id, strain_key):
    for c in ["MIC", "mic"]:
        if c in df.columns: return c
    pats = {"ecoli": ["E. coli","ecoli","e_coli"], "saureus": ["S. aureus","saureus","s_aureus","staph"]}
    for col in df.columns:
        for p in pats.get(strain_key, []):
            if p.lower() in col.lower(): return col
    return None

def load_mic_scatter(model_id, strain_key, truth_path):
    ip = INFERENCE_DIR / model_id / "battleamp-all" / "predictions.tsv"
    if not ip.exists() or not truth_path.exists(): return None
    truth = _normalise_seq_col(pd.read_csv(truth_path)); pred = _normalise_seq_col(pd.read_csv(ip, sep="\t"))
    if "sequence" not in truth.columns or "sequence" not in pred.columns: return None
    tmc = "MIC"
    if tmc not in truth.columns:
        for c in truth.columns:
            if c.lower() == "mic": tmc = c; break
    pmc = _find_mic_col(pred, model_id, strain_key)
    if pmc is None: return None
    ts = truth[["sequence", tmc]].rename(columns={tmc: "true_mic"})
    ts["true_mic"] = pd.to_numeric(ts["true_mic"], errors="coerce")
    ps = pred[["sequence", pmc]].rename(columns={pmc: "pred_mic"})
    ps["pred_mic"] = pd.to_numeric(ps["pred_mic"], errors="coerce")
    mg = ts.merge(ps, on="sequence").dropna()
    mg = mg[(mg["true_mic"] > 0) & (mg["pred_mic"] > 0)]
    return mg[["true_mic", "pred_mic"]] if len(mg) > 0 else None

def compute_mic_cache():
    MIC_CACHE.mkdir(parents=True, exist_ok=True)
    for mid in REG_FIG4:
        for sk, si in STRAIN_TRUTH.items():
            df = load_mic_scatter(mid, sk, si["path"])
            out = MIC_CACHE / f"mic_{mid}_{sk}.tsv"
            if df is not None: df.to_csv(out, sep="\t", index=False)
            else: pd.DataFrame(columns=["true_mic","pred_mic"]).to_csv(out, sep="\t", index=False)
    log.info("MIC scatter cache complete.")

_MIC_CACHE_FILES = [MIC_CACHE / f"mic_{mid}_{sk}.tsv"
                    for mid in REG_FIG4 for sk in STRAIN_TRUTH]
if not all(f.exists() for f in _MIC_CACHE_FILES):
    log.info("MIC scatter cache missing, computing...")
    compute_mic_cache()

for strain_key, strain_info in STRAIN_TRUTH.items():
    n_mod = len(REG_FIG4); n_cols = 4; n_rows = math.ceil(n_mod / n_cols)
    fig_sc, axes_sc = plt.subplots(n_rows, n_cols, figsize=(MAX_WIDTH, MAX_HEIGHT),
                                    dpi=DPI, squeeze=False)
    fig_sc.subplots_adjust(hspace=0.50, wspace=0.38, left=0.08, right=0.97, top=0.93, bottom=0.05)
    fig_sc.suptitle(f"True vs predicted MIC: {strain_info['label']}",
                    fontsize=9, fontweight="bold", y=0.97)

    mids = list(REG_FIG4.keys()); mnames = list(REG_FIG4.values())
    ticks = [1, 4, 16, 64, 256]; lims = [0.8, 600]

    for idx in range(n_rows * n_cols):
        r, c = divmod(idx, n_cols); ax = axes_sc[r, c]
        if idx >= n_mod: ax.set_visible(False); continue
        mid = mids[idx]; ax.set_facecolor(BG)
        ax.grid(True, color=GRID_CLR, lw=0.3, zorder=0); ax.set_axisbelow(True)

        cp = MIC_CACHE / f"mic_{mid}_{strain_key}.tsv"
        has = False
        if cp.exists():
            df = pd.read_csv(cp, sep="\t")
            if len(df) > 0: has = True
        if has:
            tm = df["true_mic"].values; pm = df["pred_mic"].values
            rho, _ = spearmanr(np.log2(tm), np.log2(pm))
            ax.scatter(tm, pm, s=2, alpha=0.35, color="0.30", ec="none", rasterized=True, zorder=3)
            ax.text(0.05, 0.95, f"$\\rho$ = {rho:.2f}\nn = {len(df):,}",
                    transform=ax.transAxes, ha="left", va="top", fontsize=FS_TICK,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=0.85, lw=0.4))
        else:
            ax.text(0.5, 0.5, "NO DATA", transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="red", alpha=0.4)

        ax.plot(lims, lims, "--", color="#CC2936", lw=0.5, zorder=5)
        ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xticks(ticks); ax.set_yticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks]); ax.set_yticklabels([str(t) for t in ticks])
        ax.set_title(mnames[idx], fontsize=FS_TITLE, pad=3)
        if c == 0: ax.set_ylabel("Predicted MIC", fontsize=FS_LABEL)
        if r == n_rows - 1 or idx + n_cols >= n_mod:
            ax.set_xlabel("True MIC", fontsize=FS_LABEL)

    save_figure(fig_sc, f"figure_supp_mic_{strain_info['suffix']}")
    plt.show()
