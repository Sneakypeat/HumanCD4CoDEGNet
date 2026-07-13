"""Rebuild HumanCD4CoDEGNet.cx2 WITH an embedded layout + per-element styling, so it opens
on NDEx already laid-out and coloured (and still fully interactive/draggable)."""
import pandas as pd, numpy as np, networkx as nx, os, json
from ndex2.cx2 import CX2Network, RawCX2NetworkFactory
NET=os.path.expanduser("~/CoDEG_Tcell/network")
nodes=pd.read_csv(f"{NET}/nodes.csv"); edges=pd.read_csv(f"{NET}/edges.csv")

G=nx.DiGraph()
for _,r in nodes.iterrows(): G.add_node(r["gene"],**r.to_dict())
for _,r in edges.iterrows(): G.add_edge(r["source"],r["target"],**r.to_dict())
pos=nx.spring_layout(G,k=1.7,iterations=200,seed=42)   # SAME layout as the preview PNG

COL={"gained":"#E8802E","stable":"#2E7FB0","lost":"#7A8699","peripheral":"#A9B2C0"}
od=nodes.set_index("gene")["outdeg_Stim48hr"].astype(float); odmax=od.max()
def nsize(g): return float(24+46*(od[g]/odmax))

cx=CX2Network()
idmap={}
for _,r in nodes.iterrows():
    g=r["gene"]; x,y=pos[g]
    attrs={k:(bool(r[k]) if str(r[k]) in ("True","False") else
              (int(r[k]) if isinstance(r[k],(np.integer,)) or (isinstance(r[k],float) and float(r[k]).is_integer()) else r[k]))
           for k in ["name","hub_status","is_disease_gene","is_disease_candidate","role",
                     "outdeg_Rest","outdeg_Stim8hr","outdeg_Stim48hr",
                     "indeg_Rest","indeg_Stim8hr","indeg_Stim48hr"] if k in r}
    nid=cx.add_node(attributes=attrs, x=float(x*260), y=float(-y*260))
    idmap[g]=nid
    sz=nsize(g); dis=bool(r["is_disease_candidate"])
    cx.add_node_bypass(nid, {
        "NODE_BACKGROUND_COLOR":COL[r["hub_status"]], "NODE_FILL_COLOR":COL[r["hub_status"]],
        "NODE_WIDTH":sz,"NODE_HEIGHT":sz,"NODE_SIZE":sz,
        "NODE_LABEL":g,"NODE_LABEL_COLOR":"#1B1F26","NODE_LABEL_FONT_SIZE":9,
        "NODE_BORDER_WIDTH":(5.0 if dis else 0.5),
        "NODE_BORDER_COLOR":("#E23B3B" if dis else "#5A6472"),
        "NODE_BORDER_PAINT":("#E23B3B" if dis else "#5A6472")})
for _,r in edges.iterrows():
    g=bool(r["gained_on_activation"]); col=("#E8802E" if g else "#B4BBC6")
    eid=cx.add_edge(source=idmap[r["source"]], target=idmap[r["target"]],
                    attributes={"interaction":"regulates","states":r["states"],
                                "n_states":int(r["n_states"]),"gained_on_activation":g})
    cx.add_edge_bypass(eid, {"EDGE_LINE_COLOR":col,"EDGE_STROKE_UNSELECTED_PAINT":col,
                             "EDGE_WIDTH":1.0,"EDGE_OPACITY":(85 if g else 45)})

# best-effort default style (guarded)
try:
    cx.set_visual_properties([{"default":{
        "network":{"NETWORK_BACKGROUND_COLOR":"#FFFFFF"},
        "node":{"NODE_SHAPE":"ellipse","NODE_BACKGROUND_COLOR":"#A9B2C0","NODE_LABEL_COLOR":"#1B1F26",
                "NODE_LABEL_FONT_SIZE":9,"NODE_WIDTH":30,"NODE_HEIGHT":30,"NODE_BORDER_WIDTH":0.5},
        "edge":{"EDGE_WIDTH":1.0,"EDGE_LINE_COLOR":"#B4BBC6","EDGE_OPACITY":45}},
        "nodeMapping":{},"edgeMapping":{}}])
except Exception as e: print("  (default-style note:",e,")")

NAME="HumanCD4CoDEGNet — trans-regulatory hub core of the human CD4+ T cell"
DESC=("Causal trans-regulatory hub-core network of primary human CD4+ T cells from the Marson-lab "
      "genome-scale CRISPRi Perturb-seq atlas (~22M cells; Zhu, Dann et al. 2025). Node colour = hub_status "
      "(orange gained / blue stable / grey lost on activation); size = out-degree; red border = druggable "
      "disease candidate. Edge A->B = knockdown of A significantly shifts B; orange = gained on activation. "
      "Poster https://sneakypeat.github.io/HumanCD4CoDEGNet/  Code https://github.com/Sneakypeat/HumanCD4CoDEGNet")
for k,v in {"name":NAME,"description":DESC,"version":"1.1","author":"Syed Sabih ur Rehman",
            "networkType":"causal trans-regulatory network"}.items(): cx.add_network_attribute(k,v)

out=f"{NET}/HumanCD4CoDEGNet.cx2"; cx.write_as_raw_cx2(out)
# validate
raw=open(out).read(); rt=RawCX2NetworkFactory().get_cx2network(out)
print(json.dumps({"nodes":len(rt.get_nodes()),"edges":len(rt.get_edges()),
                  "has_coords":('"x"' in raw and '"y"' in raw),
                  "has_nodeBypasses":("Bypass" in raw or "bypass" in raw),
                  "has_visualProperties":("visualProperties" in raw),
                  "bytes":len(raw)},indent=2))
print("wrote",out)
