"""Re-export GraphML carrying layout coords + color/size columns (for Cytoscape passthrough mapping)."""
import pandas as pd, numpy as np, networkx as nx, os
NET=os.path.expanduser("~/CoDEG_Tcell/network")
nodes=pd.read_csv(f"{NET}/nodes.csv"); edges=pd.read_csv(f"{NET}/edges.csv")
G=nx.DiGraph()
for _,r in nodes.iterrows(): G.add_node(r["gene"],**{k:r[k] for k in nodes.columns})
for _,r in edges.iterrows(): G.add_edge(r["source"],r["target"],**{k:r[k] for k in edges.columns})
pos=nx.spring_layout(G,k=1.7,iterations=200,seed=42)
COL={"gained":"#E8802E","stable":"#2E7FB0","lost":"#7A8699","peripheral":"#A9B2C0"}
od=nodes.set_index("gene")["outdeg_Stim48hr"].astype(float); odmax=od.max()
for g in G.nodes:
    x,y=pos[g]
    G.nodes[g]["x"]=float(x*400); G.nodes[g]["y"]=float(-y*400)     # Cytoscape Y grows down
    G.nodes[g]["color"]="#E23B3B" if bool(G.nodes[g]["is_disease_candidate"]) else COL[G.nodes[g]["hub_status"]]
    G.nodes[g]["fill"]=COL[G.nodes[g]["hub_status"]]
    G.nodes[g]["size"]=float(24+46*(od[g]/odmax))
for u,v in G.edges:
    G.edges[u,v]["color"]="#E8802E" if bool(G.edges[u,v]["gained_on_activation"]) else "#B4BBC6"
# graphml can't hold python bool objects cleanly -> cast bools to str
for n in G.nodes:
    for k,val in list(G.nodes[n].items()):
        if isinstance(val,(np.bool_,bool)): G.nodes[n][k]=str(bool(val))
        elif isinstance(val,(np.integer,)): G.nodes[n][k]=int(val)
        elif isinstance(val,(np.floating,)): G.nodes[n][k]=float(val)
for u,v in G.edges:
    for k,val in list(G.edges[u,v].items()):
        if isinstance(val,(np.bool_,bool)): G.edges[u,v][k]=str(bool(val))
        elif isinstance(val,(np.integer,)): G.edges[u,v][k]=int(val)
        elif isinstance(val,(np.floating,)): G.edges[u,v][k]=float(val)
out=f"{NET}/HumanCD4CoDEGNet.graphml"; nx.write_graphml(G,out,named_key_ids=True)
print("wrote",out,"| nodes",G.number_of_nodes(),"edges",G.number_of_edges(),
      "| now carries: x,y,color,fill,size")
