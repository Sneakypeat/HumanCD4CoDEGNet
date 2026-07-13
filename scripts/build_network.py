"""Build the HumanCD4CoDEGNet trans-regulatory HUB-CORE network and export CX2/GraphML for NDEx.

Nodes  = top out-degree regulators (validated KD) per state  ∪  disease candidates.
Edges  = causal edge regulator->target (KD of source significantly shifts target),
         reconstructed from artifacts/B_masked.npz among the node set, tagged by the
         activation state(s) in which the edge is present (Rest / Stim8hr / Stim48hr).
"""
import numpy as np, scipy.sparse as sp, pandas as pd, h5py, fsspec, json, os, time
import networkx as nx

ART=os.path.expanduser("~/CoDEG_Tcell/artifacts")
OUT=os.path.expanduser("~/CoDEG_Tcell/network"); os.makedirs(OUT,exist_ok=True)
URL="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad"
CONDS=["Rest","Stim8hr","Stim48hr"]
DISEASE=["ZAP70","ITK","LCK","PTPRC","IL12RB2"]

def dec(a): return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def obs_col(o,k):
    x=o[k]
    if isinstance(x,h5py.Group): return dec(x["categories"][:])[x["codes"][:]]
    return x[:]

# ---- load B + label maps ----
d=np.load(f"{ART}/B_masked.npz",allow_pickle=True); ncol=int(d["shape"][1])
B=np.unpackbits(d["B"],axis=1)[:,:ncol].astype(bool); keep_idx=d["keep_idx"]
gene_name=d["gene_name"].astype(str); gpos={g:i for i,g in enumerate(gene_name)}

def stream_obs():
    for attempt in range(5):
        try:
            h=h5py.File(fsspec.open(URL,block_size=4*1024*1024).open(),"r"); o=h["obs"]
            return obs_col(o,"culture_condition"), obs_col(o,"target_contrast_gene_name").astype(str)
        except Exception as e:
            if attempt==4: raise
            print(f"  stream attempt {attempt+1} failed ({type(e).__name__}); retrying...",flush=True)
            time.sleep(6*(attempt+1))
print("streaming .obs from public S3 (no bulk download) ...",flush=True)
cond_all,targ_all=stream_obs()
row_cond=cond_all[keep_idx]; row_targ=targ_all[keep_idx]
print(f"  obs mapped: {len(row_targ)} perturbation-rows",flush=True)

# ---- node selection from precomputed out-degree resource (validated KD) ----
arch=pd.read_csv(f"{ART}/arch_perturbation_outdegree.csv")   # target,gene,outdeg_*,onsig_*
arch=arch.dropna(subset=["gene"]).drop_duplicates("gene")
def top_hubs(state,k):
    s=arch[arch[f"onsig_{state}"]==True].sort_values(f"outdeg_{state}",ascending=False).head(k)
    return list(s["gene"])
hubs_rest=set(top_hubs("Rest",30)); hubs_s8=set(top_hubs("Stim8hr",30)); hubs_s48=set(top_hubs("Stim48hr",50))
nodeset=sorted((hubs_rest|hubs_s8|hubs_s48|set(DISEASE)) & set(gpos))   # must be measurable targets
print(f"  node set: {len(nodeset)} genes (Rest{len(hubs_rest)} S8{len(hubs_s8)} S48{len(hubs_s48)} +disease)",flush=True)

# ---- reconstruct edges among node set, per state (only node-set regulators -> fast) ----
nsset=set(nodeset)
edge_states={}                     # (src,dst) -> set(states)
for c in CONDS:
    m=(row_cond==c) & np.isin(row_targ,nodeset)
    subB=B[m]; subT=row_targ[m]
    for g in np.unique(subT):
        vec=subB[subT==g].max(0)   # union significant targets for regulator g in state c
        for dst in nodeset:
            if dst==g: continue
            if vec[gpos[dst]]:
                edge_states.setdefault((g,dst),set()).add(c)
    print(f"  {c}: cumulative edges {len(edge_states)}",flush=True)

# ---- node attributes ----
pg=pd.read_csv(f"{ART}/per_gene_full.csv").drop_duplicates("gene_name").set_index("gene_name")
arch_i=arch.set_index("gene")
def geti(df,idx,col,default=0):
    try: return int(round(float(df.loc[idx,col])))
    except Exception: return default
# disease gene set (best-effort from results json) + candidates
disease_genes=set(DISEASE)
try:
    dj=json.load(open(f"{ART}/disease_hubs_results.json"))
    for v in dj.values():
        if isinstance(v,list): disease_genes|={str(x) for x in v if isinstance(x,str) and x.isupper()}
        if isinstance(v,dict):
            for vv in v.values():
                if isinstance(vv,list): disease_genes|={str(x) for x in vv if isinstance(x,str) and x.isupper()}
except Exception as e: print("  (disease json note:",e,")")

def hub_status(g):
    r,s=g in hubs_rest,g in hubs_s48
    return "stable" if (r and s) else "gained" if s else "lost" if r else "peripheral"

G=nx.DiGraph()
for g in nodeset:
    o48=geti(arch_i,g,"outdeg_Stim48hr"); o0=geti(arch_i,g,"outdeg_Rest")
    i48=geti(pg,g,"nreg_Stim48hr"); i0=geti(pg,g,"nreg_Rest")
    G.add_node(g, name=g,
               outdeg_Rest=o0, outdeg_Stim8hr=geti(arch_i,g,"outdeg_Stim8hr"), outdeg_Stim48hr=o48,
               indeg_Rest=i0, indeg_Stim8hr=geti(pg,g,"nreg_Stim8hr"), indeg_Stim48hr=i48,
               hub_status=hub_status(g),
               is_disease_gene=bool(g in disease_genes),
               is_disease_candidate=bool(g in DISEASE),
               role=("broadcaster" if o48>=i48 else "receiver"))
for (s,t),st in edge_states.items():
    G.add_edge(s,t, interaction="regulates",
               state_Rest=bool("Rest" in st), state_Stim8hr=bool("Stim8hr" in st),
               state_Stim48hr=bool("Stim48hr" in st), n_states=int(len(st)),
               states=";".join([c for c in CONDS if c in st]),
               gained_on_activation=bool("Stim48hr" in st and "Rest" not in st))
# drop isolated nodes (no edges within core) so the network is clean
iso=[n for n in G.nodes if G.degree(n)==0]; G.remove_nodes_from(iso)
print(f"  final graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} directed edges "
      f"({len(iso)} isolated dropped)",flush=True)

NAME="HumanCD4CoDEGNet — trans-regulatory hub core of the human CD4+ T cell"
DESC=("Causal trans-regulatory hub-core network of primary human CD4+ T cells, reconstructed from the "
      "Marson-lab genome-scale CRISPRi Perturb-seq atlas (~22M cells; Zhu, Dann et al. 2025). Nodes are the "
      "top out-degree regulators (validated knockdown) plus druggable disease-gene candidates; a directed edge "
      "A->B means knockdown of A significantly shifts B. Edges are annotated with the activation state(s) "
      "(Rest / Stim8hr / Stim48hr) in which they are present, capturing the shape-invariant but identity-labile "
      "rewiring (hub-dominance pinned at Gini~0.92 while up to 59% of hubs turn over). "
      "Method lineage: YeastCoDEGNet (Nasar, Rehman, Ott & Alam, NAR 2026). "
      "Poster: https://sneakypeat.github.io/HumanCD4CoDEGNet/  Code: https://github.com/Sneakypeat/HumanCD4CoDEGNet")

# ---- exports ----
nx.write_graphml(G, f"{OUT}/HumanCD4CoDEGNet.graphml")
pd.DataFrame([{"gene":n,**G.nodes[n]} for n in G.nodes]).to_csv(f"{OUT}/nodes.csv",index=False)
pd.DataFrame([{"source":s,"target":t,**G.edges[s,t]} for s,t in G.edges]).to_csv(f"{OUT}/edges.csv",index=False)

import ndex2
# CX2 (modern NDEx format)
try:
    from ndex2.cx2 import NetworkXToCX2NetworkFactory
    cx2=NetworkXToCX2NetworkFactory().get_cx2network(G)
    for k,v in {"name":NAME,"description":DESC,"version":"1.0","author":"Syed Sabih ur Rehman",
                "reference":"Zhu, Dann et al. 2025 bioRxiv 10.64898/2025.12.23.696273; YeastCoDEGNet NAR 2026 gkaf1410",
                "networkType":"causal trans-regulatory network"}.items():
        cx2.add_network_attribute(k,v)
    cx2.write_as_raw_cx2(f"{OUT}/HumanCD4CoDEGNet.cx2")
    print("  wrote CX2:",f"{OUT}/HumanCD4CoDEGNet.cx2",flush=True)
except Exception as e:
    print("  CX2 export note:",repr(e),flush=True)
# CX1 (NiceCX) fallback — also NDEx-uploadable
try:
    nice=ndex2.create_nice_cx_from_networkx(G)
    nice.set_name(NAME); nice.set_network_attribute("description",DESC)
    nice.set_network_attribute("author","Syed Sabih ur Rehman")
    json.dump(nice.to_cx(), open(f"{OUT}/HumanCD4CoDEGNet.cx","w"))
    print("  wrote CX1:",f"{OUT}/HumanCD4CoDEGNet.cx",flush=True)
except Exception as e:
    print("  CX1 export note:",repr(e),flush=True)

# summary
byst={s:0 for s in ["stable","gained","lost","peripheral"]}
for n in G.nodes: byst[G.nodes[n]["hub_status"]]+=1
gained_edges=sum(1 for s,t in G.edges if G.edges[s,t]["gained_on_activation"])
print(json.dumps({"nodes":G.number_of_nodes(),"edges":G.number_of_edges(),
                  "hub_status":byst,"edges_gained_on_activation":gained_edges,
                  "disease_nodes":sorted([n for n in G.nodes if G.nodes[n]["is_disease_candidate"]])},indent=2))
print("DONE")
