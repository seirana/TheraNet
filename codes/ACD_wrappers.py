import numpy as np
import random
import read_write as rw
import networkx as nx
import pandas as pd

FINITE_INFINITY = 20000

def matrix_to_network_guney(file):

    db = rw.read_csv(file)

    G = nx.Graph()

    genes = db.loc[:,'gene1']
    genes = genes.drop_duplicates()
    genes = genes.sort_values()
    genes.index = range(len(genes))    

    for i in range(len(genes)):
        G.add_node(genes[i])

    for i in range(len(db)):
        G.add_edge(db.loc[i,'gene1'], db.loc[i,'gene2'])
        
    return G   



def matrix_to_network_edited(file, pc_genes):

    G = nx.Graph()

    pcg = rw.read_csv(pc_genes)
    db = rw.read_csv(file)

    genes = pd.concat([pcg.loc[:,'Gene stable ID'],db.loc[:,'gene1']])
    genes = genes.drop_duplicates()
    genes = genes.sort_values()
    genes.index = range(len(genes))

    for i in range(len(genes)):
        G.add_node(genes[i])

    for i in range(len(db)):
        G.add_edge(db.loc[i,'gene1'], db.loc[i,'gene2'])        
    return G    
  

  
def _bidirectional_pred_succ(G, source, target):

    if target == source:
        return ({target: None}, {source: None}, source)

    Gpred = G.adj
    Gsucc = G.adj

    # predecesssor and successors in search
    pred = {source: None}
    succ = {target: None}

    # initialize fringes, start with forward
    forward_fringe = [source]
    reverse_fringe = [target]
    while forward_fringe and reverse_fringe:
        if len(forward_fringe) <= len(reverse_fringe):
            this_level = forward_fringe
            forward_fringe = []
            for v in this_level:
                for w in Gsucc[v]:
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w] = v
                    if w in succ:  # path found
                        return pred, succ, w
        else:
            this_level = reverse_fringe
            reverse_fringe = []
            for v in this_level:
                for w in Gpred[v]:
                    if w not in succ:
                        succ[w] = v
                        reverse_fringe.append(w)
                    if w in pred:  # found path
                        return pred, succ, w
    return(FINITE_INFINITY)



def bidirectional_shortest_path(G, source, target):

    if source not in G or target not in G:
        msg = f"Either source {source} or target {target} is not in G"
        raise nx.NodeNotFound(msg)

    # call helper to do the real work
    results = _bidirectional_pred_succ(G, source, target)
    if results == FINITE_INFINITY:
        return FINITE_INFINITY
    
    pred, succ, w = results
    
    # build path from pred+w+succ
    path = []
    # from source to w
    while w is not None:
        path.append(w)
        w = pred[w]
    path.reverse()
    # from w to target
    w = succ[path[-1]]
    while w is not None:
        path.append(w)
        w = succ[w]
    return path



def shortest_path_length(G, source, target):
  p = bidirectional_shortest_path(G, source, target)
  if p == FINITE_INFINITY:
      return FINITE_INFINITY
  else:
      paths = len(p) - 1
  return paths



def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.iteritems(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = degree_to_nodes.keys()
    values = sorted(values)
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print i, low, high, len(val) 
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins



def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes



def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                #nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20): # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [ random.choice(nodes) ]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(sorted(nodes), len(nodes_selected))
        values.append(nodes_random)
    return values



def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size) 
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random



def calculate_closest_distance(network, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            #val = nx.shortest_path_length(network, node_from, node_to)
            val = shortest_path_length(network, node_from, node_to)
            values.append(val)
        d = min(values)
        values_outer.append(d)
    return values_outer



def calculate_proximity(network, nodes_from, nodes_to, n_random):
    min_bin_size = 100
    seed = 452456
    lengths = None

    nodes_network = network.nodes()

    nodes_from_pruned = set(nodes_from) & set(nodes_network)  #union
    nodes_to_pruned = set(nodes_to) & set(nodes_network)  # union
    
    if len(nodes_from_pruned) == 0 or len(nodes_to_pruned) == 0:
        return FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY, len(nodes_from_pruned), len(nodes_to_pruned), FINITE_INFINITY
    
    shortest_distances = calculate_closest_distance(network, nodes_from_pruned, nodes_to_pruned)

    all_distances = shortest_distances
    d = np.mean(shortest_distances)

    all_distances = [list(all_distances)]

    bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    nodes_from_random = get_random_nodes(nodes_from_pruned, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)

    nodes_to_random = []
    ntr = list(nodes_to_pruned)
    for i in range(n_random):
        nodes_to_random.append(ntr)

    values = np.empty(n_random)
    for i in range(n_random):
        to = []
        if type(nodes_to_random[i]) == int:
            to.append(nodes_to_random[i])
        else:
            to = nodes_to_random[i]

        frm = []
        if type(nodes_from_random[i]) == int:
            frm.append(nodes_from_random[i])
        else:
            frm = nodes_from_random[i]

        shortest_distances = calculate_closest_distance(network, frm, to)
        values[i] = np.mean(shortest_distances)

    pval = float(sum(values <= d)) / len(values)
    m, s = np.mean(values), np.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, pval, m, s, len(nodes_from_pruned), len(nodes_to_pruned), all_distances



def calculate_proximity_multiple(network, drugs_lst_file, drugs_file, disease_file, check_file, folder, sampling, title):

    drug_to_targets = rw.read_csv(drugs_file)
    drugs = rw.read_csv(drugs_lst_file)
    disease_to_genes = rw.read_csv(disease_file)
    check =  rw.read_csv(check_file)

    output = pd.DataFrame(index = range(len(drug_to_targets)),columns = ['drug','shortest distance','z score','p-value','mean','SD','number of drug target genes with PPI','number of disease target genes with PPI', 'distances', 'number of drug target genes','number of disease target genes'])

    for i in range(len(drug_to_targets)):       
        ind = [j for j, x in enumerate(drug_to_targets.loc[i,:]) if x > 0]
        nodes_from = list(drug_to_targets.iloc[i,ind])
        drug = drugs.loc[i,'drug_name']
        nodes_to = list(disease_to_genes.loc[:,'ENSEMBL ID'])
        
        if int(check.iloc[i,0]) == i:
            d, z, pval, m, s, len_nodes_from, len_nodes_to, all_distances = calculate_proximity(network, nodes_from, nodes_to, sampling)
            output.loc[i,:] = [drug,d, z, pval, m, s, len_nodes_from, len_nodes_to, all_distances, len(nodes_from), len(nodes_to)]
        else:
            output.loc[i,:] = [drug,'', '', '', '', '', '', '', '', '', '']

    file =  folder + str(sampling) + '_' + title
    rw.write_csv(output,file)