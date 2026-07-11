"""
Disease relevance of the identity-turnover: are activation-STATE-SPECIFIC ("gained") hubs enriched for
monogenic-disease / autoimmune / druggable genes vs STABLE hubs? Holds hubness constant (both are Stim-8h
hubs, both validated in Rest & 8h), guards the inducibility confound with a degree+expression-matched null
and a logistic control. Novel vs the atlas paper (which did cluster DOWNSTREAM-target GWAS enrichment; this
is regulator hubs stratified by the stable-vs-labile turnover axis).

Inputs (committed): artifacts/arch_perturbation_outdegree.csv, artifacts/per_gene_full.csv.
Gene lists fetched from the authors' public repo. Writes artifacts/disease_hubs_results.json +
artifacts/figures/fig_arch5_disease.png.
"""
import io, json, math, os, requests, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
try:
    import statsmodels.api as sm; HAVE_SM=True
except Exception: HAVE_SM=False
ART=os.path.expanduser("~/CoDEG_Tcell/artifacts"); FIG=f"{ART}/figures"; os.makedirs(FIG,exist_ok=True)
GWT="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master"
rng=np.random.default_rng(20260711)
def fetch(url): r=requests.get(url,timeout=60); r.raise_for_status(); return r.text
def listset(txt):
    out=set()
    for i,l in enumerate(txt.splitlines()):
        t=l.split("\t")[0].strip()
        if t and not (i==0 and t.lower() in {"gene","symbol","x",""}): out.add(t)
    return out
gl=lambda n: listset(fetch(f"{GWT}/metadata/gene_lists/{n}.tsv"))
ai=pd.read_csv(io.StringIO(fetch(f"{GWT}/src/6_functional_interaction/results/autoimmune_genes_opentargets.csv")))
AI={r.gene:r.max_score for r in ai.itertuples()}; AIHIGH={k for k,v in AI.items() if v>=0.5}
CLINVAR=gl("clinvar_path_likelypath"); GWAS=gl("gwascatalog")
DRUG=set().union(*[gl(c) for c in ["kinases","gpcr_union","ion_channels","nuclear_receptors","transporters","catalytic_receptors","cytokine_receptors"]])
SETS={"autoimmune(OT≥0.5)":AIHIGH,"ClinVar-path":CLINVAR,"GWAS-catalog":GWAS,"druggable":DRUG}

hub=pd.read_csv(f"{ART}/arch_perturbation_outdegree.csv")
g=pd.read_csv(f"{ART}/per_gene_full.csv")[["gene_name","baseMean"]]
hub=hub.merge(g,left_on="gene",right_on="gene_name",how="left")
both=hub[hub.onsig_Rest.fillna(False)&hub.onsig_Stim8hr.fillna(False)].dropna(subset=["outdeg_Rest","outdeg_Stim8hr","baseMean"]).copy()
for n,S in SETS.items(): both[n]=both.gene.isin(S)
K=300
topR=set(both.nlargest(K,"outdeg_Rest").gene); top8=set(both.nlargest(K,"outdeg_Stim8hr").gene)
stable,gained,lost=topR&top8,top8-topR,topR-top8
both["od_d"]=pd.qcut(both.outdeg_Stim8hr.rank(method="first"),8,labels=False)
both["bm_d"]=pd.qcut(both.baseMean.rank(method="first"),8,labels=False)
print(f"both-validated regulators: {len(both)} | stable={len(stable)} gained={len(gained)} lost={len(lost)}")

def frac(genes,col): s=both[both.gene.isin(genes)]; return s[col].mean(),int(s[col].sum()),len(s)
def matched_null(col,nperm=3000):
    gd=both[both.gene.isin(gained)]; strata=gd.groupby(["od_d","bm_d"]).size()
    pool={k:both[(both.od_d==k[0])&(both.bm_d==k[1])] for k in strata.index}
    obs=gd[col].mean(); null=np.array([pd.concat([pool[k].sample(n,replace=True,random_state=int(rng.integers(1e9)))
        for k,n in strata.items()])[col].mean() for _ in range(nperm)])
    return obs,float(null.mean()),float((obs-null.mean())/null.std()),float(((null>=obs).sum()+1)/(nperm+1))

OUT={"universe":len(both),"n_gained":len(gained),"n_stable":len(stable)}
print("\n set                    bg%  stable% gained% lost% | Fisher OR(p)   | matched z(p)  | logit OR(p)")
for n,S in SETS.items():
    bg=both[n].mean(); s_=frac(stable,n); g_=frac(gained,n); l_=frac(lost,n)
    OR,pf=stats.fisher_exact([[g_[1],g_[2]-g_[1]],[s_[1],s_[2]-s_[1]]])
    obs,nm,z,pn=matched_null(n)
    lor=lp=float("nan")
    if HAVE_SM:
        b=both.copy(); b["gained"]=b.gene.isin(gained).astype(int); b["l10bm"]=np.log10(b.baseMean+1)
        try:
            m=sm.Logit(b[n].astype(int),sm.add_constant(b[["gained","l10bm","outdeg_Stim8hr"]])).fit(disp=0)
            lor,lp=float(np.exp(m.params["gained"])),float(m.pvalues["gained"])
        except Exception: pass
    print(f" {n:20} {bg*100:4.1f} {s_[0]*100:6.1f} {g_[0]*100:6.1f} {l_[0]*100:5.1f} | {OR:4.2f} ({pf:.2f})  | {z:+4.1f} ({pn:.2f}) | {lor:4.2f} ({lp:.2f})")
    OUT[n]=dict(bg=float(bg),stable=s_[0],gained=g_[0],lost=l_[0],fisher_OR=float(OR),fisher_p=float(pf),matched_z=z,matched_p=pn,logit_OR=lor,logit_p=lp)

c=both[both.gene.isin(gained)&both["druggable"]&(both["autoimmune(OT≥0.5)"]|both["ClinVar-path"])].copy()
c["ai_score"]=c.gene.map(AI).fillna(0); c=c.sort_values("outdeg_Stim8hr",ascending=False)
OUT["candidates"]=[dict(gene=r.gene,od_Rest=int(r.outdeg_Rest),od_Stim8hr=int(r.outdeg_Stim8hr),
    ai_score=float(r.ai_score),clinvar=bool(r.gene in CLINVAR)) for r in c.itertuples()]
print("\ncandidate state-specific control points:", ", ".join(x["gene"] for x in OUT["candidates"]))
json.dump(OUT,open(f"{ART}/disease_hubs_results.json","w"),indent=2,default=float)

# ---------------- figure ----------------
plt.rcParams.update({"savefig.dpi":300,"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"axes.titlesize":11,"axes.titleweight":"bold","legend.frameon":False})
fig,ax=plt.subplots(1,2,figsize=(11.5,4.5),gridspec_kw={"width_ratios":[1.15,1]})
order=["ClinVar-path","autoimmune(OT≥0.5)","druggable","GWAS-catalog"]; labs=["ClinVar\n(monogenic)","autoimmune\nGWAS","druggable\ngenome","GWAS\ncatalog"]
xp=np.arange(4); w=.38
ax[0].bar(xp-w/2,[OUT[s]["stable"]*100 for s in order],w,color="#8AA0AE",label="stable hubs")
ax[0].bar(xp+w/2,[OUT[s]["gained"]*100 for s in order],w,color="#C85A28",label="activation-gained hubs")
for i,s in enumerate(order):
    p=OUT[s]["fisher_p"]; star="**" if p<0.05 else ("*" if p<0.1 else "n.s.")
    ax[0].text(i,max(OUT[s]["stable"],OUT[s]["gained"])*100+1.5,f"OR {OUT[s]['fisher_OR']:.1f}\n{star}",ha="center",
               fontsize=8,color="#B23B3F" if p<0.05 else "#777",fontweight="bold" if p<0.05 else "normal")
ax[0].set_xticks(xp); ax[0].set_xticklabels(labs,fontsize=8.5); ax[0].set_ylabel("% of hubs in gene set"); ax[0].set_ylim(0,48)
ax[0].legend(loc="upper center",fontsize=8.5); ax[0].set_title("A  Disease risk concentrates in the\nstate-specific (not stable) hubs")
cand=OUT["candidates"]; logy=[math.log10(x["od_Stim8hr"]+1) for x in cand]
for i in range(1,len(logy)):
    if logy[i-1]-logy[i]<0.11: logy[i]=logy[i-1]-0.11
for j,x in enumerate(cand):
    ax[1].plot([0,1],[x["od_Rest"]+1,x["od_Stim8hr"]+1],"-o",color="#C85A28",ms=5,lw=1.6,alpha=.9)
    ax[1].plot([1,1.02],[x["od_Stim8hr"]+1,10**logy[j]],color="#ccc",lw=.6)
    ax[1].text(1.04,10**logy[j],x["gene"],va="center",fontsize=9,fontweight="bold")
ax[1].set_yscale("log"); ax[1].set_xlim(-.15,1.55); ax[1].set_xticks([0,1]); ax[1].set_xticklabels(["Rest","Stim 8h"])
ax[1].set_ylabel("out-degree + 1 (targets controlled)"); ax[1].set_title("B  State-specific control points\n(gained hub ∩ druggable ∩ disease)")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_arch5_disease.png",bbox_inches="tight"); plt.close(fig)
print("wrote disease_hubs_results.json + fig_arch5_disease.png")
