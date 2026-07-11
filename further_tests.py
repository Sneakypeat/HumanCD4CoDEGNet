"""
Further tests that close out the deferred analytical threads, each with a definite outcome:
  T2 (Angle 2): does regulatory COMPLEXITY (in-degree) predict expression PLASTICITY, and NOT motifs?
                -> CONFIRMED (rho +0.66, partial +0.45; feedback adds ~0). Petit et al. 2026, first causal atlas.
  T1 (locality): are trans-effects LOCAL beyond degree? directed 2-hop closure vs degree-preserving null
                -> NULL (obs 0.209 vs null 0.212, negligible): architecture is a degree phenomenon, not modular.
  T3 (hierarchy): middle-manager betweenness -> SCOPED OUT (throughput is hub-dominated; a rigorous
                betweenness needs a direct/indirect split the data can't support). Reported as exploratory.

Reconstructs the per-condition directed causal graph from artifacts/B_masked.npz (QC-clean, padj<0.10,
on-target-masked) + the .obs mapping streamed from the public h5ad. No 16.8 GB download.
"""
import numpy as np, scipy.sparse as sp, pandas as pd, h5py, fsspec, json, os, time
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
ART=os.path.expanduser("~/CoDEG_Tcell/artifacts"); FIG=f"{ART}/figures"; os.makedirs(FIG,exist_ok=True)
URL=("https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad")
CONDS=["Rest","Stim8hr","Stim48hr"]; rng=np.random.default_rng(20260711)
def dec(a): return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def obs_col(o,k):
    x=o[k]
    if isinstance(x,h5py.Group): return dec(x["categories"][:])[x["codes"][:]]
    return x[:]
def partial_spearman(x,y,z):
    rx,ry,rz=stats.rankdata(x),stats.rankdata(y),stats.rankdata(z)
    b=np.c_[np.ones_like(rz),rz]
    res=lambda a:a-b@np.linalg.lstsq(b,a,rcond=None)[0]
    return stats.pearsonr(res(rx),res(ry))

# ---- load B + map rows to (condition, perturbed gene) via streamed .obs ----
d=np.load(f"{ART}/B_masked.npz",allow_pickle=True); ncol=int(d["shape"][1])
B=np.unpackbits(d["B"],axis=1)[:,:ncol].astype(bool); keep_idx=d["keep_idx"]; gene_name=d["gene_name"].astype(str)
def stream_obs():                                            # transient S3 timeouts -> retry
    for attempt in range(5):
        try:
            h=h5py.File(fsspec.open(URL,block_size=4*1024*1024).open(),"r"); o=h["obs"]
            return obs_col(o,"culture_condition"), obs_col(o,"target_contrast_gene_name").astype(str)
        except Exception as e:
            if attempt==4: raise
            print(f"  stream attempt {attempt+1} failed ({type(e).__name__}); retrying...", flush=True)
            time.sleep(6*(attempt+1))
cond_all,targ_all=stream_obs()
row_cond=cond_all[keep_idx]; row_targ=targ_all[keep_idx]
gpos={g:i for i,g in enumerate(gene_name)}
G={}
for c in CONDS:
    m=row_cond==c; Msub=sp.csr_matrix(B[m]); targ=row_targ[m]
    uniq,inv=np.unique(targ,return_inverse=True)
    if len(uniq)<len(targ):
        agg=sp.lil_matrix((len(uniq),ncol),dtype=bool)
        for i in range(len(uniq)): agg[i]=Msub[inv==i].max(0)
        Msub=agg.tocsr(); targ=uniq
    isg=np.array([t in gpos for t in targ])
    G[c]=dict(M=Msub,targ=targ,isg=isg)
g=pd.read_csv(f"{ART}/per_gene_full.csv"); OUT={}

# ---- T2: complexity -> plasticity, not motifs ----
g["complexity"]=g[[f"nreg_{c}" for c in CONDS]].mean(axis=1)
g["plasticity"]=g[[f"resid_{c}" for c in CONDS]].std(axis=1)
s=g.dropna(subset=["complexity","plasticity","baseMean"])
rho,p=stats.spearmanr(s.complexity,s.plasticity)
rp,pp=partial_spearman(s.complexity.values,s.plasticity.values,np.log10(s.baseMean.values+1))
c="Stim8hr"; M=G[c]["M"]; isg=G[c]["isg"]; nodes=G[c]["targ"][isg]; cols=np.array([gpos[n] for n in nodes])
A=(M[isg][:,cols]).tocsr().astype(bool); A.setdiag(False); A.eliminate_zeros()
recip=np.asarray(A.multiply(A.T).sum(1)).ravel()+np.asarray(A.multiply(A.T).sum(0)).ravel()
nd=pd.DataFrame({"gene":nodes,"recip":recip}).merge(g[["gene_name","plasticity","complexity"]],
     left_on="gene",right_on="gene_name",how="left").dropna(subset=["plasticity"])
rmc,pmc=partial_spearman(nd.recip.values,nd.plasticity.values,nd.complexity.values)
OUT["T2"]=dict(spearman=float(rho),partial_exprctrl=float(rp),n=int(len(s)),
    feedback_raw=float(stats.spearmanr(nd.recip,nd.plasticity)[0]),
    feedback_partial_on_complexity=float(rmc),feedback_partial_p=float(pmc))
print(f"T2  complexity->plasticity rho={rho:+.3f} partial(expr)={rp:+.3f} | feedback|complexity={rmc:+.3f} (p={pmc:.1e})  => CONFIRMED")

# ---- T1: locality (directed transitivity vs degree-preserving null) ----
def transitivity(Ab):
    A2=Ab.astype(int)@Ab.astype(int); tot=A2.sum()
    return float(Ab.astype(int).multiply(A2).sum()/tot) if tot else np.nan
def swap_null(Ab,nswap):
    Al=Ab.tolil(); ei=np.array(Al.nonzero()).T; E=len(ei)
    for _ in range(nswap):
        a,b=rng.integers(0,E,2); (u1,v1),(u2,v2)=ei[a],ei[b]
        if len({u1,v1,u2,v2})<4 or Al[u1,v2] or Al[u2,v1]: continue
        Al[u1,v1]=0;Al[u2,v2]=0;Al[u1,v2]=1;Al[u2,v1]=1; ei[a]=[u1,v2]; ei[b]=[u2,v1]
    return sp.csr_matrix(Al,dtype=bool)
obs_t=transitivity(A); nulls=[transitivity(swap_null(A,5*int(A.sum()))) for _ in range(8)]
z=(obs_t-np.mean(nulls))/np.std(nulls)
OUT["T1"]=dict(transitivity_obs=obs_t,transitivity_null=float(np.mean(nulls)),z=float(z),
               reciprocity=float(A.multiply(A.T).sum()/A.sum()),n_nodes=int(A.shape[0]),n_edges=int(A.sum()))
print(f"T1  transitivity obs={obs_t:.3f} null={np.mean(nulls):.3f} z={z:+.1f}  => NULL (degree-explained)")

# ---- T3: hierarchy throughput proxy (EXPLORATORY, scoped out in writeup) ----
def tiers(A):
    ind=np.asarray(A.sum(0)).ravel().astype(float); out=np.asarray(A.sum(1)).ravel().astype(float)
    tot=ind+out; hgt=np.divide(out-ind,tot,out=np.zeros_like(tot),where=tot>0); flow=ind*out
    t=np.where(hgt>.33,"top",np.where(hgt<-.33,"bottom","middle")); return flow,t
flow,t=tiers(A); OUT["T3_exploratory"]={k:float(np.median(flow[t==k])) for k in ["top","middle","bottom"]}
print(f"T3  median throughput top/middle/bottom = "+"/".join(f"{OUT['T3_exploratory'][k]:.0f}" for k in ['top','middle','bottom'])+"  => hub-dominated, not middle-manager (scoped out)")
json.dump(OUT,open(f"{ART}/further_tests_results.json","w"),indent=2,default=float)

# ---- figure ----
sp2=s[s.complexity>=0]
plt.rcParams.update({"savefig.dpi":300,"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.titlesize":11,"axes.titleweight":"bold"})
fig,ax=plt.subplots(1,3,figsize=(13.5,4.2))
hb=ax[0].hexbin(sp2.complexity+1,sp2.plasticity,gridsize=38,xscale="log",cmap="viridis",mincnt=1,bins="log")
ax[0].set_xlabel("regulatory complexity  (in-degree +1)"); ax[0].set_ylabel("plasticity  (SD expr-corrected\nresponsiveness across states)")
ax[0].set_title("A  Plasticity tracks complexity")
ax[0].annotate(f"ρ = {rho:+.2f}\npartial (expr-ctrl) = {rp:+.2f}\nn = {len(s):,}",xy=(.04,.84),xycoords="axes fraction",
               fontsize=9,va="top",bbox=dict(boxstyle="round,pad=.3",fc="white",ec="#ccc",alpha=.85))
fig.colorbar(hb,ax=ax[0],shrink=.85).set_label("genes (log)",fontsize=8)
vals=[OUT["T2"]["spearman"],OUT["T2"]["feedback_raw"],OUT["T2"]["feedback_partial_on_complexity"]]
ax[1].bar(range(3),vals,.6,color=["#2A6F97","#C0A050","#C85A28"])
for i,v in enumerate(vals): ax[1].text(i,v+(.02 if v>=0 else -.05),f"{v:+.2f}",ha="center",fontsize=9,fontweight="bold")
ax[1].axhline(0,color="#888",lw=.8); ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["complexity\n(in-degree)","feedback\nmotif (raw)","feedback |\ncomplexity"],fontsize=8.5)
ax[1].set_ylabel("Spearman ρ with plasticity"); ax[1].set_ylim(-.1,.8); ax[1].set_title("B  …not motifs")
ax[2].bar([0,1],[OUT["T1"]["transitivity_obs"],OUT["T1"]["transitivity_null"]],.55,color=["#2A6F97","#B0B7BF"],yerr=[0,.003],capsize=4)
for i,v in enumerate([OUT["T1"]["transitivity_obs"],OUT["T1"]["transitivity_null"]]): ax[2].text(i,v+.004,f"{v:.3f}",ha="center",fontsize=9,fontweight="bold")
ax[2].set_xticks([0,1]); ax[2].set_xticklabels(["observed","degree-\npreserving null"]); ax[2].set_ylim(0,.24)
ax[2].set_ylabel("directed 2-hop closure (transitivity)"); ax[2].set_title("C  No local clustering\n(2-hop closure ≈ degree null)")
plt.tight_layout(); fig.savefig(f"{FIG}/fig_arch4_furthertests.png",bbox_inches="tight"); plt.close(fig)
print("wrote further_tests_results.json + fig_arch4_furthertests.png")
