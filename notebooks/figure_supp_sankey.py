"""Supplementary: Sankey dataset overlap."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *
from matplotlib.path import Path as MplPath

def load_labels(path):
    d={}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            d[row["sequence"].strip()] = row.get("label","").strip()
    return d

amp_l = load_labels(TASKS_DIR/"amp"/"labels.tsv")
gen_l = load_labels(TASKS_DIR/"broad_activity"/"labels.tsv")
amp_pos_s = {s for s,l in amp_l.items() if l=="AMP"}
gen_pos_s = {s for s,l in gen_l.items() if l=="AMP"}; gen_neg_s = {s for s,l in gen_l.items() if l=="non-AMP"}
dbaasp_only = set(gen_l.keys()) - set(amp_l.keys())
n_amp_pos = sum(1 for l in amp_l.values() if l=="AMP"); n_amp_neg = sum(1 for l in amp_l.values() if l=="non-AMP")
D = {"n_amp_pos":n_amp_pos,"n_amp_neg":n_amp_neg,"n_gen_act":len(gen_pos_s),"n_gen_inact":len(gen_neg_s),
     "n_dbaasp_only":len(dbaasp_only),"amp_pos_to_active":len(amp_pos_s&gen_pos_s),
     "amp_pos_to_inactive":len(amp_pos_s&gen_neg_s),"dbaasp_only_active":len(dbaasp_only&gen_pos_s),
     "dbaasp_only_inactive":len(dbaasp_only&gen_neg_s)}

C_AP="#3A76AF";C_AN="#9DC3E6";C_GA=ACTIVE_COLOR;C_GI=INACTIVE_COLOR;C_DB="#8C7BAE";C_LC="#C44E52"
def _bz(ax,x0,y0t,y0b,x1,y1t,y1b,clr,z=1):
    cx=(x0+x1)/2
    v=[(x0,y0t),(cx,y0t),(cx,y1t),(x1,y1t),(x1,y1b),(cx,y1b),(cx,y0b),(x0,y0b),(x0,y0t)]
    c=[MplPath.MOVETO]+[MplPath.CURVE4]*3+[MplPath.LINETO]+[MplPath.CURVE4]*3+[MplPath.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(MplPath(v,c),fc=clr,ec="none",lw=0,zorder=z))

fsk,axk=plt.subplots(figsize=(MAX_WIDTH,2.5), dpi=DPI); axk.set_xlim(-0.01,1.01); axk.set_ylim(-0.05,1.05); axk.axis("off")
nw=0.022; xL=0.30; xR=0.68; ref=max(D["n_amp_pos"],D["n_amp_neg"])
sc=0.22/(ref**0.35); gap=0.06
def nh(n): return max(sc*(n**0.35),0.025) if n>0 else 0
h_ap=nh(D["n_amp_pos"]);h_an=nh(D["n_amp_neg"]);h_db=nh(D["n_dbaasp_only"]);h_ga=nh(D["n_gen_act"]);h_gi=nh(D["n_gen_inact"])
af=[D["amp_pos_to_active"],D["amp_pos_to_inactive"],D["dbaasp_only_active"],D["dbaasp_only_inactive"]]
sf=0.08/max(af)
def fh(n): return max(sf*n,0.003)

yt=0.92
y_ap_t=yt;y_ap_b=yt-h_ap;y_an_t=y_ap_b-gap;y_an_b=y_an_t-h_an;y_db_t=y_an_b-gap;y_db_b=y_db_t-h_db
y_ga_t=yt;y_ga_b=yt-h_ga;y_gi_t=y_ga_b-gap;y_gi_b=y_gi_t-h_gi
def nd(x,yb,w,ht,c):
    axk.add_patch(mpatches.FancyBboxPatch((x,yb),w,ht,boxstyle="round,pad=0.002",fc=c,ec="white",lw=0.4,zorder=5))
nd(xL,y_ap_b,nw,h_ap,C_AP);nd(xL,y_an_b,nw,h_an,C_AN);nd(xL,y_db_b,nw,h_db,C_DB)
nd(xR,y_ga_b,nw,h_ga,C_GA);nd(xR,y_gi_b,nw,h_gi,C_GI)

pad=0.01
for t,b,nm,n,c in [(y_ap_t,y_ap_b,"AMP",D["n_amp_pos"],"black"),(y_an_t,y_an_b,"non-AMP",D["n_amp_neg"],"black"),
    (y_db_t,y_db_b,"DBAASP-only",D["n_dbaasp_only"],C_DB)]:
    axk.text(xL-pad,(t+b)/2,f"{nm}  (n={n:,})",fontsize=FS_TICK,ha="right",va="center",fontweight="bold",color=c)
for t,b,nm,n in [(y_ga_t,y_ga_b,"Active",D["n_gen_act"]),(y_gi_t,y_gi_b,"Inactive",D["n_gen_inact"])]:
    axk.text(xR+nw+pad,(t+b)/2,f"{nm}  (n={n:,})",fontsize=FS_TICK,ha="left",va="center",fontweight="bold")

sp=0.001; al=0.30
fhc=fh(D["amp_pos_to_active"]);fhx=fh(D["amp_pos_to_inactive"]);fhda=fh(D["dbaasp_only_active"]);fhdi=fh(D["dbaasp_only_inactive"])
s1t=y_ap_t;s1b=s1t-fhc;d1t=y_ga_t;d1b=d1t-fhc
_bz(axk,xL+nw,s1t,s1b,xR,d1t,d1b,(*mcolors.to_rgb(C_GA),al),2)
s2t=s1b-sp;s2b=s2t-fhx;d2t=y_gi_t;d2b=d2t-fhx
_bz(axk,xL+nw,s2t,s2b,xR,d2t,d2b,(*mcolors.to_rgb(C_LC),0.45),3)
s3t=y_db_t;s3b=s3t-fhda;d3t=d1b-sp;d3b=d3t-fhda
_bz(axk,xL+nw,s3t,s3b,xR,d3t,d3b,(*mcolors.to_rgb(C_DB),al),2)
s4t=s3b-sp;s4b=s4t-fhdi;d4t=d2b-sp;d4b=d4t-fhdi
_bz(axk,xL+nw,s4t,s4b,xR,d4t,d4b,(*mcolors.to_rgb(C_DB),al),2)

mx=(xL+nw+xR)/2
def fm(*a): return sum(a)/len(a)
axk.text(mx,fm(s1t,s1b,d1t,d1b),f"{D['amp_pos_to_active']:,} (consistent)",fontsize=FS_TICK,ha="center",va="center",color="#2B7A3E",fontstyle="italic",zorder=10)
axk.text(mx,fm(s2t,s2b,d2t,d2b),f"{D['amp_pos_to_inactive']:,} (label change)",fontsize=FS_TICK,ha="center",va="center",color="#8B1A1A",fontweight="bold",zorder=10)
axk.text(mx+0.015,fm(s3t,s3b,d3t,d3b),f"{D['dbaasp_only_active']:,}",fontsize=FS_TICK,ha="center",va="center",color="#6B5F8A",zorder=10)
axk.text(mx+0.015,fm(s4t,s4b,d4t,d4b),f"{D['dbaasp_only_inactive']:,}",fontsize=FS_TICK,ha="center",va="center",color="#6B5F8A",zorder=10)
axk.text(xL+nw/2,0.99,"AMP/non-AMP dataset",fontsize=FS_LABEL,ha="center",va="bottom",fontweight="bold",color=C_AP)
axk.text(xR+nw/2,0.99,"GeneralActivity dataset",fontsize=FS_LABEL,ha="center",va="bottom",fontweight="bold",color="#444")

axk.legend(handles=[mpatches.Patch(fc=C_AP,label="AMP"),mpatches.Patch(fc=C_GA,label=r"Active (MIC $\leq$ 32 $\mu$g/ml)"),
    mpatches.Patch(fc=C_AN,label="non-AMP"),mpatches.Patch(fc=C_GI,label=r"Inactive (MIC $\geq$ 128 $\mu$g/ml)"),
    mpatches.Patch(fc=C_DB,label="DBAASP-only"),mpatches.Patch(fc=mcolors.to_rgba(C_LC,0.45),ec=C_LC,lw=0.5,label="Label change")],
    loc="lower center",bbox_to_anchor=(0.48,-0.03),ncol=3,fontsize=FS_TICK,frameon=False,
    handlelength=1,handleheight=0.7,columnspacing=0.8,labelspacing=0.25)
fsk.tight_layout(pad=0.1); save_figure(fsk,"sankey_dataset_overlap"); plt.show()
