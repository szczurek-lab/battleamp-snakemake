"""Figure 5: MD simulations and ML predictor comparison."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, pearsonr

safe5 = pd.read_csv(MD_DIR / "safe" / "summary_aa_dual.csv")
paired5 = pd.read_csv(MD_DIR / "paired" / "summary_aa_dual.csv")
safe_avg = safe5.groupby(["sequence","is_active","metric","step"])["value"].mean().reset_index()
paired_avg = paired5.groupby(["sequence","is_active","metric","step"])["value"].mean().reset_index()

safe_ml = pd.read_csv(RESULTS_ROOT / "aggregated" / "control-md_predictions.tsv", sep="\t")
y_safe5 = (safe_ml["true_class"] == "AMP").values

paired_gt = pd.read_csv(MD_DIR / "activity-cliffs.tsv", sep="\t")
paired_gt["MIC_clean"] = paired_gt["MIC [ug/ml]"].astype(str).str.replace(",","").astype(float)
paired_pred = pd.read_csv(RESULTS_ROOT / "aggregated" / "activity-cliffs_predictions.tsv", sep="\t")
pair_labels = {}
for i in range(len(paired_gt)//2):
    pair_labels[paired_gt.iloc[2*i]["sequence"].upper()] = True
    pair_labels[paired_gt.iloc[2*i+1]["sequence"].upper()] = False

STEPS = [1,2,3,4,5]

def plot_split_violins(ax, df, mc="S", yl="S-score"):
    data = df[df["metric"]==mc].copy(); data["Status"]=data["is_active"].map({1:"Active",0:"Inactive"})
    steps = sorted(data["step"].unique())
    for st, off, clr in [("Active",-0.17,ACTIVE_COLOR),("Inactive",0.17,INACTIVE_COLOR)]:
        sub = data[data["Status"]==st]
        dbs = [sub[sub["step"]==s]["value"].dropna().values for s in steps]
        if all(len(d)==0 for d in dbs): continue
        vp = ax.violinplot(dbs, positions=np.array(steps)+off, widths=0.3, showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(clr); body.set_alpha(0.85); body.set_edgecolor("none"); body.set_linewidth(0)
            m = np.mean(body.get_paths()[0].vertices[:,0])
            if off<0: body.get_paths()[0].vertices[:,0] = np.clip(body.get_paths()[0].vertices[:,0], -np.inf, m)
            else: body.get_paths()[0].vertices[:,0] = np.clip(body.get_paths()[0].vertices[:,0], m, np.inf)
        vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(0.5)
    for y in np.arange(0,1.1,0.2): ax.axhline(y, color=GRID_CLR, lw=0.3, zorder=0)
    ax.set_xticks(steps); ax.set_xlabel("Simulation block (100 ns each)")
    ax.set_ylabel(yl); ax.set_xlim(0.4, 5.6); style_ax(ax)

def compute_auc_mat(df, metrics, steps):
    mat = np.full((len(steps), len(metrics)), np.nan)
    for j, mc in enumerate(metrics):
        for i, st in enumerate(steps):
            sub = df[(df["metric"]==mc)&(df["step"]==st)]
            if len(sub)<5 or len(np.unique(sub["is_active"]))<2: continue
            try: mat[i,j] = roc_auc_score(sub["is_active"].values, sub["value"].values)
            except: pass
    for j in range(len(metrics)):
        c=mat[:,j]; v=c[~np.isnan(c)]
        if len(v)>0 and np.mean(v)<0.5: mat[:,j]=1-c
    return mat

def plot_auc_hm(ax, mat, metrics, steps, vmin=0.35, vmax=0.90, ml=None,
                bur_sep=None, hel_sep=None, hel_end=None):
    im = ax.imshow(mat, cmap="RdYlGn", aspect="auto", vmin=vmin, vmax=vmax, interpolation="nearest")
    for i in range(len(steps)):
        for j in range(len(metrics)):
            v=mat[i,j]
            if not np.isnan(v):
                ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=FS_TICK,color="white" if v>0.82 or v<0.42 else "#444")
    ax.set_xticks(range(len(metrics))); ax.set_xticklabels(ml or metrics, fontsize=FS_TICK)
    ax.set_yticks(range(len(steps))); ax.set_yticklabels([f"Block {s}" for s in steps])
    ax.set_ylabel("AUC")
    nm = len(metrics)
    y_bar, y_tick, y_lbl = -0.26, -0.18, -0.34
    def _bracket(x0, x1, label):
        ax.plot([x0, x1], [y_bar, y_bar], transform=ax.transAxes, color="#555", lw=0.7, clip_on=False)
        for xb in [x0, x1]:
            ax.plot([xb, xb], [y_tick, y_bar], transform=ax.transAxes, color="#555", lw=0.7, clip_on=False)
        ax.text((x0+x1)/2, y_lbl, label, ha="center", va="top", fontsize=FS_TICK,
                transform=ax.transAxes, clip_on=False)
    if hel_sep is not None:
        ax.axvline(hel_sep - 0.5, color="white", lw=1.5, zorder=3)
        he = hel_end if hel_end is not None else (bur_sep or nm)
        _bracket(hel_sep / nm, he / nm, "Helicity")
    if bur_sep is not None:
        ax.axvline(bur_sep - 0.5, color="white", lw=1.5, zorder=3)
        _bracket(bur_sep / nm, 1.0, "BUR")
    return im

def plot_zp(ax, df, step, legend=True):
    zn=[f"ZP-{b}" for b in ["0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5"]]
    zp=[f"ZP{b}" for b in ["0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5"]]
    def ga(bins):
        a=[]
        for mc in bins:
            sub=df[(df["metric"]==mc)&(df["step"]==step)]
            if len(sub)<5 or len(np.unique(sub["is_active"]))<2: a.append(np.nan); continue
            try: a.append(roc_auc_score(sub["is_active"].values, sub["value"].values))
            except: a.append(np.nan)
        return a
    an=ga(zn); ap=ga(zp); x=np.arange(len(zn)); w=0.35
    ax.bar(x-w/2,an,w,color=ZP_NEG_COLOR,alpha=0.75,ec="none",label="ZP neg. (thinning)")
    ax.bar(x+w/2,ap,w,color=ZP_POS_COLOR,alpha=0.75,ec="none",label="ZP pos. (thickening)")
    ax.axhline(0.5,color="black",lw=0.6,ls="--",alpha=0.6,label="AUC = 0.5")
    for y in [0.3,0.4,0.6,0.7,0.8]: ax.axhline(y,color=GRID_CLR,lw=0.3,zorder=0)
    ax.set_xticks(x); ax.set_xticklabels([f"{v:.1f}" for v in [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]], fontsize=FS_TICK)
    ax.set_xlabel(f"Phosphate displacement (Å), block {step}"); ax.set_ylabel("AUC"); ax.set_ylim(0.15,0.92)
    if legend:
        ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="0.7", fontsize=FS_LEGEND)
    style_ax(ax)

def plot_ml_box(ax, dd, models):
    n=len(models); pa=np.arange(n)-0.18; pi=np.arange(n)+0.18; bw=0.30
    np.random.seed(42)
    for idx,(mk,_,_,dn) in enumerate(models):
        act,ina=dd.get(mk,([],[]))
        if len(act)==0 or len(ina)==0: continue
        al=np.log2(np.clip(act,0.01,None)); il=np.log2(np.clip(ina,0.01,None))
        for v,p,c in [(al,pa[idx],ACTIVE_COLOR),(il,pi[idx],INACTIVE_COLOR)]:
            ax.boxplot([v],positions=[p],widths=bw,patch_artist=True,showfliers=False,zorder=2,
                boxprops=dict(fc=c,lw=0,ec="none",alpha=0.9),medianprops=dict(color="black",lw=0.5),
                whiskerprops=dict(color="black",lw=0.5),capprops=dict(color="black",lw=0.5))
            ax.scatter(p+np.random.normal(0,0.04,len(v)),v,s=6,color=c,ec="black",lw=0.2,alpha=0.8,zorder=3)
    ax.set_xticks(np.arange(n)); ax.set_xticklabels([m[3] for m in models],fontsize=FS_TICK,rotation=35,ha="right")
    ax.set_xlim(-0.6,n-0.4); ax.set_ylabel(r"$\log_2$ MIC (ML regressor predictions)"); style_ax(ax)

ml_safe = {}
for mk,col,_,_ in ML_MODELS_FIG5:
    mask=safe_ml[col].notna(); vals=safe_ml.loc[mask,col].values; lab=y_safe5[mask.values]
    ml_safe[mk] = (vals[lab], vals[~lab])
ml_cliff = {}
for mk,_,cv,_ in ML_MODELS_FIG5:
    sub=paired_pred[(paired_pred["variant"]==cv)&paired_pred["MIC"].notna()]; am,im=[],[]
    for _,r in sub.iterrows():
        seq=str(r["sequence"]).upper()
        if seq in pair_labels:
            if pair_labels[seq]: am.append(r["MIC"])
            else: im.append(r["MIC"])
    ml_cliff[mk] = (np.array(am), np.array(im))

rho_safe_d = {}
for mk,col,_,_ in ML_MODELS_FIG5:
    mask=safe_ml[col].notna()
    if not mask.any(): continue
    v=safe_ml.loc[mask,col].values; lb=y_safe5[mask.values].astype(int)
    rho_safe_d[mk] = spearmanr(v, lb)[0]

rho_cliff_d = {}
for mk,_,cv,_ in ML_MODELS_FIG5:
    sub=paired_pred[(paired_pred["variant"]==cv)&paired_pred["MIC"].notna()]
    vv,ll=[],[]
    for _,r in sub.iterrows():
        seq=str(r["sequence"]).upper()
        if seq in pair_labels:
            vv.append(r["MIC"]); ll.append(1 if pair_labels[seq] else 0)
    if len(vv)>=5: rho_cliff_d[mk]=spearmanr(vv,ll)[0]

fig5 = plt.figure(figsize=(MAX_WIDTH, MAX_HEIGHT), dpi=DPI)
fig5.text(0.27, 0.97, "Control", ha="center", va="top", fontsize=FS_TITLE, fontweight="bold")
fig5.text(0.78, 0.97, "Activity cliffs", ha="center", va="top", fontsize=FS_TITLE, fontweight="bold")
gs5_top = gridspec.GridSpec(3, 2, hspace=0.55, wspace=0.35, left=0.09, right=0.96, top=0.93, bottom=0.32,
                            height_ratios=[1, 1.15, 1.05])
gs5_bot = gridspec.GridSpec(1, 2, wspace=0.35, left=0.09, right=0.96, top=0.20, bottom=0.06)

LMD = [Patch(fc=ACTIVE_COLOR,alpha=0.85,ec="none",label="Active"),
       Patch(fc=INACTIVE_COLOR,alpha=0.85,ec="none",label="Inactive")]
ax5a=fig5.add_subplot(gs5_top[0,0]); plot_split_violins(ax5a,safe_avg,"S","S-score")
ax5a.legend(handles=LMD,loc="lower right",frameon=True,facecolor="white",edgecolor="0.7")
ax5b=fig5.add_subplot(gs5_top[0,1]); plot_split_violins(ax5b,paired_avg,"S","S-score")
ax5b.legend(handles=LMD,loc="lower right",frameon=True,facecolor="white",edgecolor="0.7")

for ax_v,df_v in [(ax5a,safe_avg),(ax5b,paired_avg)]:
    for st in STEPS:
        sub=df_v[(df_v["metric"]=="S")&(df_v["step"]==st)]
        if len(sub)>=5 and len(np.unique(sub["is_active"]))>1:
            r,_=pearsonr(sub["value"].values,sub["is_active"].values.astype(float))
            ax_v.text(st,1.02,f"r={r:.2f}",ha="center",va="bottom",fontsize=FS_ANNOT,
                      transform=ax_v.get_xaxis_transform(),clip_on=False)

hm5=["S","H1","H9","HM","H5","WAT_BUR_REL","WAT_BUR_MAX","PEPT_BUR_MAX"]
hml=["S","H1","H9","HM","H5","Wat\nRel","Wat\nMax","Pept\nMax"]
sa5=compute_auc_mat(safe_avg,hm5,STEPS); pa5=compute_auc_mat(paired_avg,hm5,STEPS)
aav=np.concatenate([sa5[~np.isnan(sa5)],pa5[~np.isnan(pa5)]])
vmi=max(0.35,np.floor(aav.min()*20)/20); vma=min(0.92,np.ceil(aav.max()*20)/20)
plot_auc_hm(fig5.add_subplot(gs5_top[1,0]),sa5,hm5,STEPS,vmi,vma,hml,bur_sep=5,hel_sep=1,hel_end=5)
plot_auc_hm(fig5.add_subplot(gs5_top[1,1]),pa5,hm5,STEPS,vmi,vma,hml,bur_sep=5,hel_sep=1,hel_end=5)
plot_zp(fig5.add_subplot(gs5_top[2,0]),safe_avg,5,False)
plot_zp(fig5.add_subplot(gs5_top[2,1]),paired_avg,5)

ax5g=fig5.add_subplot(gs5_bot[0,0]); plot_ml_box(ax5g,ml_safe,ML_MODELS_FIG5)
ax5h=fig5.add_subplot(gs5_bot[0,1]); plot_ml_box(ax5h,ml_cliff,ML_MODELS_FIG5)
allml=[]
for d in [ml_safe,ml_cliff]:
    for a,i in d.values():
        if len(a)>0: allml.extend(np.log2(np.clip(a,0.01,None)))
        if len(i)>0: allml.extend(np.log2(np.clip(i,0.01,None)))
ymi=np.percentile(allml,0.5)-0.5; yma=np.percentile(allml,99.5)+0.5
ax5g.set_ylim(ymi,yma); ax5h.set_ylim(ymi,yma); ax5g.invert_yaxis(); ax5h.invert_yaxis()
ax5g.legend(handles=[Patch(fc=ACTIVE_COLOR,lw=0,alpha=0.9,label="Active"),
    Patch(fc=INACTIVE_COLOR,lw=0,alpha=0.9,label="Inactive")],
    fontsize=FS_LEGEND,loc="lower left",frameon=True,facecolor="white",edgecolor="0.7",handlelength=1,handleheight=0.6)

for ax_m,rho_d in [(ax5g,rho_safe_d),(ax5h,rho_cliff_d)]:
    tr=ax_m.get_xaxis_transform()
    for idx,(mk,_,_,_) in enumerate(ML_MODELS_FIG5):
        if mk in rho_d:
            ax_m.text(idx,1.02,f"ρ={rho_d[mk]:.2f}",ha="center",va="bottom",fontsize=FS_ANNOT,
                      transform=tr,clip_on=False,rotation=45)

for lbl,ax in [("a",ax5a),("b",fig5.axes[2]),("c",fig5.axes[4]),("d",ax5g),
               ("e",ax5b),("f",fig5.axes[3]),("g",fig5.axes[5]),("h",ax5h)]:
    p=ax.get_position()
    xoff = 0.07
    fig5.text(p.x0-xoff,p.y1+0.01,lbl,fontsize=FS_PANEL,fontweight="bold",va="top",transform=fig5.transFigure)

save_figure(fig5, "figure5")
plt.show()
