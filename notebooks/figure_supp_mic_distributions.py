"""Supplementary: MIC distributions."""
import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *

ACTIVITY_DIR = DATA_ROOT / "activity"
ACT_THR = 32; INACT_THR = 128
ACT_L2 = np.log2(ACT_THR); INACT_L2 = np.log2(INACT_THR)
L2_MIN = -7; L2_MAX = 12
BIN_EDGES = np.arange(L2_MIN, L2_MAX+1, dtype=float)
XTICK_L2 = np.array([-6,-4,-2,0,2,4,6,8,10,12], dtype=float)

DATASETS_MIC = [
    ("a","GeneralActivity","csv","broad.csv",(3415,940)),
    ("b","Gram+","csv","gramplus.csv",(2064,1136)),
    ("c","Gram$-$","csv","gramminus.csv",(2937,1075)),
    ("d",r"$\it{S.}$ $\it{aureus}$","filter",("gramplus.csv","species","Staphylococcus aureus"),(2064,1136)),
    ("e",r"$\it{E.}$ $\it{coli}$","filter",("gramminus.csv","species","Escherichia coli"),(2586,1025)),
    ("f",r"$\it{A.}$ $\it{baumannii}$","filter",("gramminus.csv","species","Acinetobacter baumannii"),(653,221)),
    ("g",r"$\it{K.}$ $\it{pneumoniae}$","filter",("gramminus.csv","species","Klebsiella pneumoniae"),(572,383)),
    ("h",r"$\it{P.}$ $\it{aeruginosa}$","filter",("gramminus.csv","species","Pseudomonas aeruginosa"),(1380,868)),
    ("i",r"$\it{S.}$ $\it{aureus}$ ATCC 25923","csv",os.path.join("strain","staphylococcusaureusatcc25923.csv"),(858,391)),
    ("j",r"$\it{S.}$ $\it{aureus}$ ATCC 33591","csv",os.path.join("strain","staphylococcusaureusatcc33591.csv"),(47,51)),
    ("k",r"$\it{S.}$ $\it{aureus}$ ATCC 43300","csv",os.path.join("strain","staphylococcusaureusatcc43300.csv"),(238,94)),
    ("l",r"$\it{E.}$ $\it{coli}$ ATCC 25922","csv",os.path.join("strain","escherichiacoliatcc25922.csv"),(1445,517)),
    ("m",r"$\it{A.}$ $\it{baumannii}$ ATCC 19606","csv",os.path.join("strain","acinetobacterbaumanniiatcc19606.csv"),(278,100)),
    ("n",r"$\it{K.}$ $\it{pneumoniae}$ ATCC 700603","csv",os.path.join("strain","klebsiellapneumoniaeatcc700603.csv"),(214,115)),
    ("o",r"$\it{P.}$ $\it{aeruginosa}$ ATCC 27853","csv",os.path.join("strain","pseudomonasaeruginosaatcc27853.csv"),(699,392)),
]

def load_min_mic(dd, st, sa):
    if st=="csv": df=pd.read_csv(os.path.join(dd,sa))
    else: cp,col,val=sa; df=pd.read_csv(os.path.join(dd,cp)); df=df[df[col]==val]
    mm=df.groupby("id")["activity"].min(); return mm.clip(lower=2.0**L2_MIN)

def bar_clrs(be):
    c=[]
    for i in range(len(be)-1):
        if be[i+1]<=ACT_L2: c.append(ACTIVE_COLOR)
        elif be[i]>=INACT_L2: c.append(INACTIVE_COLOR)
        else: c.append("#999999")
    return c

nR,nC=5,3
fmic,axmic=plt.subplots(nR,nC,figsize=(MAX_WIDTH,7.8),dpi=DPI); fmic.subplots_adjust(left=0.09,right=0.97,bottom=0.08,top=0.97,hspace=0.55,wspace=0.35)
bc=bar_clrs(BIN_EDGES); centres=0.5*(BIN_EDGES[:-1]+BIN_EDGES[1:])
dd=str(ACTIVITY_DIR)

for idx,(pl,dn,st,sa,exp) in enumerate(DATASETS_MIC):
    ax=axmic.flat[idx]; mm=load_min_mic(dd,st,sa); l2=np.log2(mm.values)
    na=int((mm<=ACT_THR).sum()); ni=int((mm>=INACT_THR).sum())
    counts,_=np.histogram(l2,bins=BIN_EDGES)
    ax.bar(centres,counts,width=0.88,color=bc,ec="white",lw=0.25,zorder=3)
    ax.axvspan(ACT_L2,INACT_L2,color=EXCLUDED_ZONE,alpha=0.5,zorder=1)
    for xp in (ACT_L2,INACT_L2): ax.axvline(xp,color="black",ls="--",lw=0.7,zorder=4)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4,integer=True))
    ax.grid(axis="y",color=GRID_CLR,lw=0.3,zorder=0); ax.set_axisbelow(True)
    ax.set_xticks(XTICK_L2); ax.set_xticklabels([f"{int(v)}" for v in XTICK_L2],fontsize=FS_TICK)
    ax.set_xlim(L2_MIN-0.5,L2_MAX+0.5); ax.set_title(f"({pl}) {dn}",fontsize=FS_LABEL,loc="left",pad=3)
    ax.text(0.03,0.95,f"n = {na+ni}\nactive: {na}\ninactive: {ni}",transform=ax.transAxes,fontsize=FS_ANNOT,
        va="top",ha="left",bbox=dict(boxstyle="round,pad=0.25",fc="white",ec="#ccc",lw=0.3,alpha=0.9),zorder=5)
    if idx%nC==0: ax.set_ylabel("Number of peptides",fontsize=FS_LABEL)
    if idx>=(nR-1)*nC: ax.set_xlabel(r"log$_2$ MIC ($\mu$g/ml)",fontsize=FS_LABEL)

fmic.legend(handles=[Patch(fc=ACTIVE_COLOR,ec="white",lw=0.3,label=r"Active (MIC $\leq$ 32 $\mu$g/ml)"),
    Patch(fc="#999999",ec="white",lw=0.3,label="Excluded intermediate"),
    Patch(fc=INACTIVE_COLOR,ec="white",lw=0.3,label=r"Inactive (MIC $\geq$ 128 $\mu$g/ml)")],
    loc="lower center",ncol=3,fontsize=FS_LEGEND,frameon=False,bbox_to_anchor=(0.5,-0.01))
save_figure(fmic,"mic_distributions"); plt.show()
