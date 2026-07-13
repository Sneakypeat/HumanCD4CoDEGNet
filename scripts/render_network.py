"""Render a preview PNG of the hub-core network (dark, on-brand)."""
import pandas as pd, numpy as np, networkx as nx, matplotlib, os
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
NET=os.path.expanduser("~/CoDEG_Tcell/network")
nodes=pd.read_csv(f"{NET}/nodes.csv"); edges=pd.read_csv(f"{NET}/edges.csv")
G=nx.DiGraph()
for _,r in nodes.iterrows(): G.add_node(r["gene"],**r.to_dict())
for _,r in edges.iterrows(): G.add_edge(r["source"],r["target"],gained=bool(r["gained_on_activation"]))
BG="#0B0E13"; COL={"gained":"#E8935A","stable":"#6FB3D6","lost":"#8A94A6","peripheral":"#454C5C"}
od=nodes.set_index("gene")["outdeg_Stim48hr"].astype(float)
sizes=[70+1000*(od[n]/od.max()) for n in G.nodes]
colors=[COL.get(G.nodes[n]["hub_status"],"#454C5C") for n in G.nodes]
pos=nx.spring_layout(G,k=1.7,iterations=200,seed=42)
fig,ax=plt.subplots(figsize=(12,9),dpi=115); fig.patch.set_facecolor(BG); ax.set_facecolor(BG); ax.axis("off")
eg=[(u,v) for u,v in G.edges if G.edges[u,v]["gained"]]; es=[(u,v) for u,v in G.edges if not G.edges[u,v]["gained"]]
nx.draw_networkx_edges(G,pos,edgelist=es,edge_color="#5A6472",width=0.5,alpha=0.10,arrows=False,ax=ax)
nx.draw_networkx_edges(G,pos,edgelist=eg,edge_color="#E8935A",width=0.6,alpha=0.13,arrows=False,ax=ax)
nx.draw_networkx_nodes(G,pos,node_size=sizes,node_color=colors,edgecolors=BG,linewidths=0.8,ax=ax)
DISEASE=["ZAP70","ITK","LCK","PTPRC","IL12RB2"]
top=list(od.sort_values(ascending=False).head(6).index)
lab={n:n for n in set(DISEASE)|set(top) if n in G.nodes}
nx.draw_networkx_labels(G,pos,labels=lab,font_size=9.5,font_color="#EAECEF",font_weight="bold",ax=ax)
# ring disease genes
nx.draw_networkx_nodes(G,pos,nodelist=[n for n in DISEASE if n in G.nodes],
    node_size=[70+1000*(od[n]/od.max())+120 for n in DISEASE if n in G.nodes],
    node_color="none",edgecolors="#ff5d5d",linewidths=1.8,ax=ax)
ax.set_title("HumanCD4CoDEGNet — trans-regulatory hub core of the human CD4$^+$ T cell",
             color="#EAECEF",fontsize=15,fontweight="bold",pad=14)
leg=[Line2D([0],[0],marker='o',color='none',markerfacecolor=COL[k],markersize=11,
      label={"gained":"hub gained on activation","stable":"stable hub","lost":"hub lost on activation","peripheral":"peripheral"}[k])
      for k in ["gained","stable","lost","peripheral"]]
leg.append(Line2D([0],[0],marker='o',color='none',markerfacecolor='none',markeredgecolor='#ff5d5d',
      markeredgewidth=1.8,markersize=13,label="druggable disease candidate"))
lg=ax.legend(handles=leg,loc="lower left",frameon=False,fontsize=9.5,labelcolor="#C7CDD6")
ax.text(0.5,-0.02,"node size ∝ out-degree (broadcast) · orange edges = gained when the cell activates",
        transform=ax.transAxes,ha="center",color="#9AA4B0",fontsize=9.5)
out=f"{NET}/HumanCD4CoDEGNet_network.png"; fig.savefig(out,facecolor=BG,dpi=115,bbox_inches="tight")
import shutil; shutil.copy(out,os.path.expanduser("~/Movies/HumanCD4CoDEGNet_network.png"))
print("wrote",out)
