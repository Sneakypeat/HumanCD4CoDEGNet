"""Rebuild a STRICTLY-VALID CX2 for NDEx: nodes (with layout coords) + edges + attributes only.
No custom visualProperties/bypasses (those tripped NDEx's strict parser). Colour in NDEx/Cytoscape
by the embedded `hub_status` / out-degree attributes (one click)."""
import pandas as pd, numpy as np, networkx as nx, os, json
from ndex2.cx2 import CX2Network, RawCX2NetworkFactory
NET=os.path.expanduser("~/CoDEG_Tcell/network")

# ---- diagnose the current (broken) file around the reported column ----
cur=f"{NET}/HumanCD4CoDEGNet.cx2"
if os.path.exists(cur):
    raw=open(cur).read()
    print("=== context around col 291206 of the OLD file ===")
    print(repr(raw[291120:291300]))
    print("=== aspects present in OLD file ===")
    try:
        arr=json.loads(raw); print([list(a.keys())[0] for a in arr if isinstance(a,dict) and a])
    except Exception as e: print("  old file json note:",e)

nodes=pd.read_csv(f"{NET}/nodes.csv"); edges=pd.read_csv(f"{NET}/edges.csv")
G=nx.DiGraph()
for _,r in nodes.iterrows(): G.add_node(r["gene"])
for _,r in edges.iterrows(): G.add_edge(r["source"],r["target"])
pos=nx.spring_layout(G,k=1.7,iterations=200,seed=42)

def clean(v):
    if isinstance(v,(np.bool_,bool)): return bool(v)
    if isinstance(v,(np.integer,)): return int(v)
    if isinstance(v,(np.floating,float)):
        f=float(v); return None if (f!=f or f in (float("inf"),float("-inf"))) else f
    return str(v)

cx=CX2Network(); idmap={}
NATTR=["name","hub_status","is_disease_gene","is_disease_candidate","role",
       "outdeg_Rest","outdeg_Stim8hr","outdeg_Stim48hr","indeg_Rest","indeg_Stim8hr","indeg_Stim48hr"]
for _,r in nodes.iterrows():
    g=r["gene"]; x,y=pos[g]
    attrs={k:clean(r[k]) for k in NATTR if k in r and clean(r[k]) is not None}
    idmap[g]=cx.add_node(attributes=attrs, x=round(float(x*300),2), y=round(float(-y*300),2))
EATTR=["interaction","states","n_states","gained_on_activation"]
for _,r in edges.iterrows():
    attrs={k:clean(r[k]) for k in EATTR if k in r and clean(r[k]) is not None}
    attrs.setdefault("interaction","regulates")
    cx.add_edge(source=idmap[r["source"]], target=idmap[r["target"]], attributes=attrs)

NAME="HumanCD4CoDEGNet — trans-regulatory hub core of the human CD4+ T cell"
DESC=("Causal trans-regulatory hub-core network of primary human CD4+ T cells from the Marson-lab "
      "genome-scale CRISPRi Perturb-seq atlas (~22M cells; Zhu, Dann et al. 2025). Nodes are the top "
      "out-degree regulators (validated knockdown) per activation state plus druggable disease-gene "
      "candidates; a directed edge A->B means knockdown of A significantly shifts B, annotated with the "
      "activation state(s) it appears in. Colour by hub_status (gained/stable/lost) and size by out-degree. "
      "Poster https://sneakypeat.github.io/HumanCD4CoDEGNet/  Code https://github.com/Sneakypeat/HumanCD4CoDEGNet")
for k,v in {"name":NAME,"description":DESC,"version":"1.2","author":"Syed Sabih ur Rehman",
            "networkType":"causal trans-regulatory network"}.items(): cx.add_network_attribute(k,v)

out=f"{NET}/HumanCD4CoDEGNet.cx2"; cx.write_as_raw_cx2(out)

# ---- STRICT validation (mirrors NDEx expectations) ----
txt=open(out).read()
bad=[]
json.loads(txt, parse_constant=lambda c: bad.append(c))   # flags NaN/Infinity/-Infinity
arr=json.loads(txt)
assert isinstance(arr,list), "top-level must be a JSON array"
aspects=[list(a.keys())[0] for a in arr if isinstance(a,dict) and a]
rt=RawCX2NetworkFactory().get_cx2network(out)
print("\n=== NEW clean file ===")
print(json.dumps({"nodes":len(rt.get_nodes()),"edges":len(rt.get_edges()),
                  "aspects":aspects,"nan_or_inf_tokens":bad,
                  "has_coords":('"x"' in txt and '"y"' in txt),
                  "has_visualProperties":("visualProperties" in txt),
                  "has_bypasses":("Bypasses" in txt),"bytes":len(txt)},indent=2))
print("STRICT-VALID" if not bad else "WARNING: NaN/Inf present")
