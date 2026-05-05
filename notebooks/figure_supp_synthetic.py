"""Supplementary: Synthetic negative properties."""
import sys
import string
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from figure_config import *
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
SYNTH_COLORS = {"Active": ACTIVE_COLOR, "Inactive": INACTIVE_COLOR,
                "Random": RANDOM_COLOR, "Realistic": REALISTIC_COLOR, "Shuffled": SHUFFLED_COLOR}
SYNTH_ORDER = ["Active","Inactive","Random","Realistic","Shuffled"]

def clean_seqs(seqs): return [s.strip() for s in seqs if all(c in VALID_AA for c in s.strip()) and len(s.strip())>=3]
def pos_charge_frac(seqs): return np.array([sum(c in "KRH" for c in s)/len(s) for s in seqs])

def comp_desc(seqs):
    seqs=clean_seqs(seqs)
    g=GlobalDescriptor(seqs);g.calculate_charge(ph=7.4);charge=g.descriptor.flatten()
    g2=GlobalDescriptor(seqs);g2.isoelectric_point();pi=g2.descriptor.flatten()
    g3=GlobalDescriptor(seqs);g3.aliphatic_index();al=g3.descriptor.flatten()
    g4=GlobalDescriptor(seqs);g4.aromaticity();ar=g4.descriptor.flatten()
    p=PeptideDescriptor(seqs,"eisenberg");p.calculate_global();hy=p.descriptor.flatten()
    p2=PeptideDescriptor(seqs,"eisenberg");p2.calculate_moment(window=1000,angle=100);hm=p2.descriptor.flatten()
    return {"length":np.array([len(s) for s in seqs],dtype=float),"charge":charge,"isoelectric_point":pi,
            "positive_charge_fraction":pos_charge_frac(seqs),"hydrophobicity":hy,"hydrophobic_moment":hm,
            "aliphatic_index":al,"aromaticity":ar}

rng=np.random.RandomState(42); NS=3000
adf=pd.read_csv(TASKS_DIR/"broad_activity"/"labels.tsv",sep="\t")
ap=adf.loc[adf["label"]=="AMP","sequence"].tolist(); an=adf.loc[adf["label"]=="non-AMP","sequence"].tolist()
syn={}
for tag in ["random","realistic","shuffled"]:
    sd=pd.read_csv(TASKS_DIR/f"synthetic_{tag}"/"labels.tsv",sep="\t")
    syn[tag.capitalize()]=sd.loc[sd["label"]=="non-AMP","sequence"].tolist()

dsets={"Active":list(rng.choice(ap,min(NS,len(ap)),replace=False)),
       "Inactive":list(rng.choice(an,min(NS,len(an)),replace=False))}
for nm in ["Random","Realistic","Shuffled"]: dsets[nm]=list(rng.choice(syn[nm],NS,replace=False))
props={nm:comp_desc(s) for nm,s in dsets.items()}

fl={"Active":adf.loc[adf["label"]=="AMP","sequence"].str.len().values.astype(float),
    "Inactive":adf.loc[adf["label"]=="non-AMP","sequence"].str.len().values.astype(float)}
for tag in ["Random","Realistic","Shuffled"]:
    sd=pd.read_csv(TASKS_DIR/f"synthetic_{tag.lower()}"/"labels.tsv",sep="\t")
    fl[tag]=sd.loc[sd["label"]=="non-AMP","sequence"].str.len().values.astype(float)

fsp,axsp=plt.subplots(2,4,figsize=(MAX_WIDTH,4),dpi=DPI)
pdefs=[("length","Length (residues)",True),("charge","Net charge (pH 7.4)",False),
    ("isoelectric_point","Isoelectric point",False),("positive_charge_fraction","Positive charge\nfraction",False),
    ("hydrophobicity","Hydrophobicity\n(Eisenberg)",False),("hydrophobic_moment","Hydrophobic moment",False),
    ("aliphatic_index","Aliphatic index",False),("aromaticity","Aromaticity",False)]

for ai,(key,yl,uf) in enumerate(pdefs):
    ax=axsp.flat[ai]; dl=[]; rl=[]
    for ds in SYNTH_ORDER:
        d=fl[ds] if (key=="length" and uf) else props[ds][key]; rl.append(d)
        lo,hi=np.percentile(d,[1,99]); dl.append(d[(d>=lo)&(d<=hi)])
    parts=ax.violinplot(dl,positions=range(len(SYNTH_ORDER)),showmeans=False,showmedians=False,showextrema=False,
        widths=0.72,bw_method="silverman",points=200)
    for body,ds in zip(parts["bodies"],SYNTH_ORDER):
        body.set_facecolor(SYNTH_COLORS[ds]);body.set_edgecolor("none");body.set_alpha(0.82)
    for i,d in enumerate(rl):
        q1,med,q3=np.percentile(d,[25,50,75])
        ax.vlines(i,q1,q3,color="k",lw=0.7); ax.scatter(i,med,color="white",s=6,zorder=3,ec="k",lw=0.4)
    ax.set_ylabel(yl,fontsize=FS_LABEL); ax.set_xticks(range(len(SYNTH_ORDER)))
    ax.set_xticklabels(SYNTH_ORDER,rotation=45,ha="right",fontsize=FS_TICK)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=5))
    style_ax(ax)
    ax.text(-0.02,1.05,string.ascii_lowercase[ai],transform=ax.transAxes,fontsize=9,fontweight="bold",va="bottom",ha="right")

fsp.legend(handles=[mpatches.Patch(fc=SYNTH_COLORS[ds],label=ds,ec="none",alpha=0.82) for ds in SYNTH_ORDER],
    loc="upper center",ncol=5,fontsize=FS_TICK,frameon=False,bbox_to_anchor=(0.5,1.01),
    handlelength=1,handletextpad=0.3,columnspacing=0.8)
plt.tight_layout(rect=[0,0,1,0.95],h_pad=1.8,w_pad=1.2)
save_figure(fsp,"fig_synthetic_properties"); plt.show()
