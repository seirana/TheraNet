#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: s.hashemi
"""

import read_write as rw
import matplotlib.pyplot as plt
import networkx as nx
import PIL
import numpy as np

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

file = '~/output/PSC_nodes_drug_list(drug_network)'
drug_lst = rw.read_csv(file)

file = '~/output/PSC_nodes_gene_list_type(drug_network)'
gene_lst_tp = rw.read_csv(file)

file = '~/output/PSC_gene_drug_gene_matrix(drug_network)'
mat = rw.read_csv(file)
mat = mat.to_numpy()

dm = mat.shape
c =[]
r =[]
m = []
for i in range(dm[1]):
    s = sum(mat[:,i])
    if s > 0:
        c.append(i)
        if i < len(drug_lst):
            m.append(i)

for i in range(dm[0]):
    s = sum(mat[i,:])
    if s > 0:
        r.append(i)

ls = np.zeros((len(drug_lst),1), dtype=int)
for i in range(len(ls)):
    ls[i] = i
 
for i in range(len(drug_lst)-1):
    if ls[i] == i:
        for j in range(i+1,len(drug_lst)):
            if ls[j] == j:
                d = 0
                for k in range(dm[0]):
                    d = d + abs(mat[k,j]-mat[k,i])
                if d == 0:
                    ls[j] = i
                    mat[:,j] = 0                  
 
# Image URLs for graph nodes
icons = {
    "drug": "~/icons/drug.png",
    "drug_and_disease_target": "~/icons/drug_and_disease_target.png",
    "drug_target": "~/icons/drug_target.png",
    "indirect_drug_target": "~/icons/indirect_drug_target.png",
    "non_drug_target": "~/icons/non_drug_target.png",
}

# Load images
images = {k: PIL.Image.open(fname) for k, fname in icons.items()}

i1 = [j for j, x in enumerate(gene_lst_tp.loc[:,'gene_type']) if x == "drug_and_disease_target"]
i2 = [j for j, x in enumerate(gene_lst_tp.loc[:,'gene_type']) if x == "drug_target"]
i3 = [j for j, x in enumerate(gene_lst_tp.loc[:,'gene_type']) if x == "indirect_drug_target"]
i4 = [j for j, x in enumerate(gene_lst_tp.loc[:,'gene_type']) if x == "non_drug_target"]

d = 0
for i in range(len(drug_lst)):
    if sum(mat[:,i]) > 0:
        d = d+1
        
subset_sizes = [d,i1,i2,i3,i4]
    
G = nx.Graph()
cnt = -1
for i in range(len(drug_lst)):
    if sum(mat[:,i]) > 0:
        label = drug_lst.loc[i,'drugs']
        G.add_node(label, image=images["drug"],layer=0)

for i in range(len(gene_lst_tp)):
    if sum(mat[i,:]) > 0:

        if isfloat(gene_lst_tp.loc[i,'gene_id']):
            label = 'nan'
            gene_lst_tp.loc[i,'gene_id'] = 'nan'
        label = gene_lst_tp.loc[i,'gene_id']

        if gene_lst_tp.loc[i,'gene_type'] == "drug_and_disease_target":
            G.add_node(label, image=images["drug_and_disease_target"],layer=1)

        if gene_lst_tp.loc[i,'gene_type'] == "drug_target":
            G.add_node(label, image=images["drug_target"], layer=2)

        if gene_lst_tp.loc[i,'gene_type'] == "indirect_drug_target":
            G.add_node(label, image=images["indirect_drug_target"], layer=3)

        if gene_lst_tp.loc[i,'gene_type'] == "non_drug_target":
            G.add_node(label, image=images["non_drug_target"], layer=4)

for i in range(len(drug_lst)):
    for j in range(len(gene_lst_tp)):
        if mat[j,i] > 0:
            labelj = gene_lst_tp.loc[j,'gene_id']
            labeli = drug_lst.loc[i,'drugs']  
            if gene_lst_tp.loc[j,'gene_type'] == 'drug_and_disease_target':      
                G.add_edge(labelj, labeli, color='pink')
            if gene_lst_tp.loc[j,'gene_type'] == 'drug_target':      
                G.add_edge(labelj, labeli, color='darkred')

for i in range(len(gene_lst_tp)):
    for j in range(i+1,len(gene_lst_tp)):
        if mat[j,i+len(drug_lst)] > 0:
            labelj = gene_lst_tp.loc[j,'gene_id']
            labeli = gene_lst_tp.loc[i,'gene_id']
            G.add_edge(labelj, labeli, color='darkblue')

pos = nx.multipartite_layout(G, subset_key="layer")
fig, ax = plt.subplots(figsize=(12,8))

node_sizes = []
for n in range(len(G.nodes)):
        node_sizes.append(1)
nx.draw_networkx_nodes(G, pos, nodelist=G.nodes, node_size=node_sizes)

edges = G.edges()
colors = [G[u][v]['color'] for u,v in edges]

nx.draw_networkx_edges(
    G,
    pos=pos,
    width=2.0, 
    alpha=0.6,
    edgelist=G.edges,
    edge_color = colors,
)

a = 0.012
for i in pos:
    pos[i] = pos[i]+a
    
nx.draw_networkx_labels(G, pos, font_size=8, font_color="black", 
                        horizontalalignment='left', 
                        verticalalignment='bottom'
)

for i in pos:
    pos[i] = pos[i]-a


# Transform from data coordinates (scaled between xlim and ylim) to display coordinates
tr_figure = ax.transData.transform
# Transform from display to figure coordinates
tr_axes = fig.transFigure.inverted().transform

# Select the size of the image (relative to the X axis)
icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.05
icon_center = icon_size / 2.0

fig.canvas.draw()

# Add the respective image to each node
for n in G.nodes:
    xf, yf = tr_figure(pos[n])
    xa, ya = tr_axes((xf, yf))
    a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
    a.imshow(G.nodes[n]["image"])
    a.axis("off")
    
plt.show()
fig.savefig('~/output/drug_network.png', dpi=1000)
