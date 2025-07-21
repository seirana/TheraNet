#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""
from bs4 import BeautifulSoup
import read_write as rw
import networkx as nx
import pandas as pd
import numpy as np
import requests
import random
import re
import os
#..........................................................................................................
#
#
#drug_networks
#
#
#..........................................................................................................   
'''
FINITE_INFINITY number must be bigger than the longest path in the network of genes. For simplicity 
assumed equal to the number of human genes (20,000)
'''
FINITE_INFINITY = 20000

def matrix_to_network_guney():
    '''
    generates gene-gene network based on Guney's algorithm
    '''
    file = './TheraNet/data/gene_gene_PPI700_ENSEMBL'
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



def matrix_to_network_TheraNet():
    '''
    generates gene-gene network based on Guney's algorithm, but the edition of 
    considering all genes, including genes with PPIs or not.
    '''
    file = './TheraNet/data/gene_gene_PPI700_ENSEMBL'
    pc_genes = './TheraNet/data/protein_coding_genes_ENSEMBL'
    pcg = rw.read_csv(pc_genes)
    db = rw.read_csv(file)

    G = nx.Graph()

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

    # predecessors and successors in search
    pred = {source: None}
    succ = {target: None}

    # Initialize fringes, start with forward
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
        if len(set(nodes_from_pruned) & set(nodes_to_pruned)) == 0:
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



def calculate_proximity_multiple(disease_file):
    
    disease_to_genes = rw.read_csv(disease_file)
    drug_to_targets = rw.read_csv('./TheraNet/data/DtoDGI_ENSEMBL(filtered)')
    drugs = rw.read_csv('./TheraNet/data/drugs(filtered)')
    check =  rw.read_csv('./TheraNet/data/repeated(filtered)')
    
    for i in range(2):
        if i == 0:
            which_method = 'TheraNet'
            network = matrix_to_network_TheraNet()
            sampling = 100000
            parts = 100
        else: 
            which_method = 'guney'
            network = matrix_to_network_guney()
            sampling = 1000
            parts = 1
        
        for part in range(parts):
            l = len(drug_to_targets)
            b = int(part*l/parts)
            e = int((part+1)*l/parts)-1
           
            output = pd.DataFrame(index = range(e-b+1),columns = ['drug','shortest distance','z score','p-value','mean','SD','number of drug target genes with PPI','number of disease target genes with PPI', 'distances', 'number of drug target genes','number of disease target genes'])
           
            for i in range(b,e+1):
                ind = [j for j, x in enumerate(drug_to_targets.loc[i,:]) if x > 0]
                nodes_from = list(drug_to_targets.iloc[i,ind])
                drug = drugs.loc[i,'drug_name']
                nodes_to = list(disease_to_genes.loc[:,'ENSEMBL ID'])
                
                if int(check.iloc[i,0]) == i:
                    d, z, pval, m, s, len_nodes_from, len_nodes_to, all_distances = calculate_proximity(network, nodes_from, nodes_to, sampling)
                    output.loc[i-b,:] = [drug,d, z, pval, m, s, len_nodes_from, len_nodes_to, all_distances, len(nodes_from), len(nodes_to)]
                else:
                    output.loc[i-b,:] = [drug,'', '', '', '', '', '', '', '', '', '']
            
            file =  './TheraNet/output/drug_scores' + which_method  + str(part) 
            if i == 1:
                './TheraNet/output/drug_scores_' + which_method
            rw.write_csv(output,file)
    #..........................................................................
    #merge files
    #..........................................................................
    sampling = 100000
    e = pd.DataFrame()
    for j in range(int(sampling/1000)):
        file = './TheraNet/output/drug_scores_TheraNet' +str(j)
        TheraNet = rw.read_csv(file)
        e = pd.concat([e,TheraNet],axis=0)
    rw.write_csv(e,'./TheraNet/output/drug_scores_TheraNet')
    
    for j in range(int(sampling/1000)):
        file = './TheraNet/output/drug_scores_TheraNet' +str(j)
        if os.path.exists(file):
            os.remove(file)

    file = './TheraNet/data/repeated(filtered)'
    rep = rw.read_csv(file)
    
    for i in range(2):
        if i == 0:
            file = './TheraNet/output/drug_scores_TheraNet'
        else:
            file = './TheraNet/output/drug_scores_guney'
            
    drugs = rw.read_csv(file)
    w = len(drugs.T)
    inf = pd.DataFrame(index=range(len(e)),columns=e.columns)
    for i in range(len(drugs)):
        if i != int(rep.iloc[i,0]):
            inf.iloc[i,1:w] = drugs.iloc[int(rep.iloc[i,0]),1:w]
            inf.iloc[i,0] = drugs.iloc[i,0]
        else:
              inf.iloc[i,:] = drugs.iloc[i,:]
    rw.write_csv(inf,file)
#..........................................................................................................
#
#
#check_for_side_effects
#
#
#..........................................................................................................       
def classify_tissue_side_effect(text, tissue_keywords):
     """
     Extracts tissue-related sentences and then searches for severity keywords within them.
     Returns one of: 'high', 'medium', 'low', 'rare', or 'none'.
     """
     # Define severity keywords.
     severity_levels = {
         'high': ['severe', 'high', 'significant', 'marked', 'atrophy of'],
         'medium': ['moderate', 'intermediate'],
         'low': ['mild', 'low', 'minimal'],
         'rare': ['rare', 'infrequent', 'seldom'],
     }
     
     problem_descrption = {
         'prefix': ['monitor', 'may cause', 'worsening','check for','fatal','drug-induced','toxic','have','abnormalities in'],
         'postfix': ['problem', 'dysfunction', 'injury','impairment','disease', 'blood tests', 'function test']
    }

     matches = []
     which = 0
     text = text.lower()
     # Patterns
     for tissue in tissue_keywords:
        tissue = tissue.lower().strip()
        
        for postfix in problem_descrption['postfix']:
            postfix = postfix.lower().strip()
            pattern = " ".join(filter(None, [tissue, postfix])).strip()
            if pattern and re.search(re.escape(pattern), text):
                matches.append(pattern)
                which = -1
                
            for prefix in problem_descrption['prefix']:
                prefix = prefix.lower().strip()
                pattern = " ".join(filter(None, [prefix, tissue, postfix])).strip()
                if pattern and re.search(re.escape(pattern), text):
                    matches.append(pattern)
                    which = 1
                    for keywords,severitys  in severity_levels.items():
                        for severity in severitys:
                            severity = severity.lower().strip()
                            pattern = " ".join(filter(None, [prefix, severity, tissue, postfix])).strip()
                            if pattern and re.search(re.escape(pattern), text):
                                matches.append(pattern)
                                which = 1
                            pattern = " ".join(filter(None, [severity, prefix, tissue, postfix])).strip()
                            if pattern and re.search(re.escape(pattern), text):
                                matches.append(pattern)
                                which = 1
     if len(matches) > 0:
         if which == -1:
             result = f'May have side-effects on the {tissue_keywords[0]}'
         if which == 1:
             result = f'It has side-effects on the {tissue_keywords[0]}'
     else:
         result = f'There is no information about side-effect on the {tissue_keywords[0]}'
     return result



def query_rx_or_drugs_com(drug_link):
    url = drug_link
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        text = soup.get_text()
        return text
    return ''



def side_effects(drug,tissue_keywords):
    drug_links = rw.read_csv('./TheraNet/data/drugs_links')
    drug_links['drug_name'] = drug_links['drug_name'].str.lower()
    drug_links = drug_links[drug_links['drug_name'].str.contains(drug, na=False)]
    side_effect = []
    if len(drug_links)>0:
        drug_link = drug_links.iloc[0,6]
        if pd.isna(drug_link):
            drug_link = f'https://www.drugs.com/mtm/{drug}.html'
    else:
        drug_link = f'https://www.drugs.com/mtm/{drug}.html'
    text = query_rx_or_drugs_com(drug_link)
    if len(text) > 0:
        side_effect = classify_tissue_side_effect(text,tissue_keywords)

    if len(side_effect) == 0:
        if len(drug_links)>0:
            drug_link = drug_links.iloc[0,5]
            if pd.isna(drug_link):
                drug_link = f'https://www.rxlist.com/{drug}-drug.htm'
        else:
            drug_link = f'https://www.rxlist.com/{drug}-drug.htm'
        text = query_rx_or_drugs_com(drug_link)
        if len(text) > 0:
            side_effect = classify_tissue_side_effect(text,tissue_keywords)
    if len(side_effect) == 0:
        side_effect = 'No information about possible side effects'        
    return side_effect



def query_drugbank_and_drug_info(drug,dosage_search):
    with open("./TheraNet/data/full database.xml", "r", encoding="utf-8") as file:
        text = file.read()

    drug_links = rw.read_csv('./TheraNet/data/drugs_links')
    drug_links['drug_name'] = drug_links['drug_name'].str.lower()
    drug_id = drug_links[drug_links['drug_name'].str.contains(drug, na=False)]
    if len(drug_id)>0:
        drug_id = drug_id.iloc[0,0]

        if pd.notna(drug_id):
            pattern = re.compile(
                rf'<drugbank-id\s+primary="true">{re.escape(drug_id)}</drugbank-id>(.*?)<drugbank-id\s+primary=',
                re.IGNORECASE | re.DOTALL
            )
            match = pattern.search(text)
            if match:
                drug_text = match.group(1).strip()
                match = re.search(r"<description>(.*?)\.", drug_text)
                if match:
                    routes = ''
                    drug_description = match.group(1).strip()
                    unique_routes  = ''
                    if dosage_search == 'True':
                        routes = re.findall(r'<route>\s*(.*?)\s*</route>', drug_text, flags=re.IGNORECASE | re.DOTALL)
                        unique_routes = list(set(routes))
                    return drug_description, unique_routes
    return '',''



def query_drugbank_and_drug_dosage(drug):
    drug_links = rw.read_csv('./TheraNet/data/drugs_links')
    drug_links['drug_name'] = drug_links['drug_name'].str.lower()
    drug_links = drug_links[drug_links['drug_name'].str.contains(drug, na=False)]
    if len(drug_links)>0:
        drug_link = drug_links.iloc[0,6]
        if pd.isna(drug_link):
            drug_link = f'https://www.drugs.com/mtm/{drug}.html'
    else:
        drug_link = f'https://www.drugs.com/mtm/{drug}.html'

    text = query_rx_or_drugs_com(drug_link)
    if len(text) > 0:
        matchs = re.search(r"Dosage form:\s(.*?)(?=\n)", text)
        if matchs:
            drug_dosage = matchs.group(1).strip()
            return drug_dosage, 'False'
    return [], 'True'



def tissue_keywords_(tissue):
    tissue = tissue.lower()
    tissue_keywords = {
        #this line must be adopted with the desired tissue name and synonyms
        'liver':['liver', 'hepatic', 'hepatotoxic', 'hepatitis','hepatotoxicity']
        }

    if len(tissue_keywords[tissue])== 0:
        return([])
    return tissue_keywords[tissue]



def aggregate_info(drug,tissue):
    tissue_keywords = tissue_keywords_(tissue)
    if len(tissue_keywords )> 0:
        drug = drug.lower()
        drug_dosage, dosage_search = query_drugbank_and_drug_dosage(drug)
        description, formulation = query_drugbank_and_drug_info(drug, dosage_search)
        side_effect = side_effects(drug,tissue_keywords)
        result = pd.DataFrame()
        result.loc[0,'drug'] = drug
        if dosage_search == 'True':
            result.loc[0,'drug dosage'] = str(formulation)
        else:
            result.loc[0,'drug dosage'] = str(drug_dosage)

        result.loc[0,'tissue_side-effect'] = str(side_effect)
        result.loc[0,'drug_description'] = description
        return result
    return ('The tissue name is not found in the database.')



def check_for_drug_info():
    drug_list = rw.read_csv('./TheraNet/output/drug_scores_TheraNet')
    tissue = 'liver' #this line must be adopted with the desired tissue name

    df = pd.DataFrame()
    for drug in drug_list['drug']:
        result = aggregate_info(drug,tissue)
        df = pd.concat([df, result], ignore_index=True)

    df = pd.concat([drug_list, df], ignore_index=True)
    df.to_csv("./TheraNet/output/drug_scores_TheraNet.csv")
#..........................................................................................................
#
#
#drug_scoring
#
#
#..........................................................................................................
def matrix_to_network(drug_genes, PPI, disease_genes, pc_genes):

    network = nx.Graph()
    all_drug_genes = pd.DataFrame(drug_genes.values.ravel(), columns=['merged'])

    genes = pd.concat([pc_genes.loc[:,'Gene stable ID'],PPI.loc[:,'gene1'],disease_genes.loc[:,'ENSEMBL ID'],all_drug_genes.iloc[:,0]])
    genes = genes.drop_duplicates()
    genes = genes.sort_values()
    genes.index = range(len(genes))

    for i in range(len(genes)):
        network.add_node(genes[i])

    for i in range(len(PPI)):
        network.add_edge(PPI.loc[i,'gene1'], PPI.loc[i,'gene2'])
    return network



def weight_a_disease_gene(network,disease_gene, PPI):
    weight = 1
    neighbors = [n for n in network.neighbors(disease_gene)]
    if len(neighbors) > 0:
        inds = [k for k, x in enumerate(PPI.loc[:,'gene1']) if x == disease_gene]
        for j in range(len(inds)):
            weight = weight + PPI.loc[inds[j],'max_ppi']
    return weight



def PPI_weight(gene1,gene2, PPI):
    inds1 = [k for k, x in enumerate(PPI.loc[:,'gene1']) if x == gene1]
    inds2 = [k for k, x in enumerate(PPI.loc[:,'gene2']) if x == gene2]
    ind = list(set(inds1) & set(inds2))
    return float(PPI.loc[ind[0],'max_ppi'])



def sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene(network,drug_gene,disease_gene, PPI):
    weight  = 0
    for path in nx.all_simple_paths(network, source=drug_gene, target=disease_gene, cutoff=2):
       weight_path = 1
       for k in range(len(path)-1):
           weight_path = weight_path * PPI_weight(path[k], path[k+1], PPI)
       weight  = weight + weight_path * weight_a_disease_gene(network,disease_gene, PPI)
    return weight



def weight_a_drug(network,drug, dgi, drug_genes, PPI, disease_genes):
    weight = 0
    clmns = [k for k, x in enumerate(drug_genes.iloc[drug,:]) if x > 0]
    for i in range(len(clmns)):
        weight_gene = 0
        for j in range(len(disease_genes)):
            weight_gene = weight_gene + sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene(network,drug_genes.iloc[drug,i],disease_genes.iloc[j,0],PPI)
        weight  = weight + weight_gene * dgi.iloc[drug,i]
    return weight



def score_drugs(network, dgi, drug_genes, drugs, PPI, disease_genes, candidate_drugs):
    weights = np.zeros((len(drugs),1), dtype=float)
    c=0
    for drug in range(len(drugs)):
        if float(candidate_drugs.loc[drug,'z score']) <= -1.96 and float(candidate_drugs.loc[drug,'p-value']) <= 0.05:
            c=c+1
            print(c)
            weights[drug,1] = weight_a_drug(network,drug, dgi, drug_genes, PPI, disease_genes)
    return weights



def drug_scoring(disease_file):

    dgi = rw.read_csv('./TheraNet/data/DtoGI_scores(filtered)')
    drug_genes = rw.read_csv('./TheraNet/data/DtoGI_ENSEMBL(filtered)')
    drugs = rw.read_csv('./TheraNet/data/drugs(filtered)')
    pc_genes = rw.read_csv('./TheraNet/data/protein_coding_genes_ENSEMBL')
    PPI = rw.read_csv('./TheraNet/data/gene_gene_PPI700_ENSEMBL')
    PPI = PPI.fillna(0)
    PPI.loc[:,'max_ppi'] = PPI.loc[:,'max_ppi']/1000

    disease_genes = rw.read_csv('./TheraNet/data/'+disease_file)
    
    file = './TheraNet/output/drug_scores_TheraNet'
    candidate_drugs = rw.read_csv(file)
    
    network = matrix_to_network(drug_genes, PPI, disease_genes, pc_genes)
    
    scores = score_drugs(network, dgi, drug_genes, drugs, PPI, disease_genes, candidate_drugs)
    score = pd.DataFrame(scores)
    score.columns = ['TheraNet_Drug_score']
    
    candidate_drugs = pd.concat([candidate_drugs, scores], ignore_index=True)
    rw.write_csv(candidate_drugs,'./TheraNet/output/scoring_drugs')
#..........................................................................................................
#
#
#disease_target_gene_scoring
#
#
#..........................................................................................................  
def sum_weight_all_paths_between_a_drug_gene_and_a_disease_gene(network,drug_gene,disease_gene, PPI):
    weight  = 0
    for path in nx.all_simple_paths(network, source=drug_gene, target=disease_gene, cutoff=2):
       weight_path = 1
       for k in range(len(path)-1):
           weight_path = weight_path * PPI_weight(path[k], path[k+1])
       weight  = weight + weight_path
    return weight



def weight_a_drug_for_disease_gene(network,drug,disease, PPI, drug_genes, disease_genes, dgi):
    weight = 0
    clmns = [k for k, x in enumerate(drug_genes.iloc[drug,:]) if x > 0]
    for i in range(len(clmns)):
        weight =  weight + sum_weight_all_paths_between_a_drug_gene_and_a_disease_gene(network,drug_genes.iloc[drug,i],disease_genes.iloc[disease,0])
    weight = weight * dgi.iloc[drug,i]
    return weight



def score_disease_genes(network,PPI, drug_genes, disease_genes, drugs, dgi, candidate_drugs): 
    weights = np.zeros((len(drugs),len(disease_genes)), dtype=float)
    c=0
    for drug in range(len(drugs)):
        if float(candidate_drugs.loc[drug,'z score']) <= -1.96 and float(candidate_drugs.loc[drug,'p-value']) <= 0.05:
            c=c+1
            print(c)
            for disease in range(len(disease_genes)): 
                weights[drug,disease] = weight_a_drug_for_disease_gene(network,drug,disease, PPI, drug_genes, disease_genes, dgi)
    return weights



def diease_target_gene_scoring(disease_file):
    dgi = rw.read_csv('./TheraNet/data/DtoGI_scores(filtered)')
    drug_genes = rw.read_csv('./TheraNet/data/DtoGI_ENSEMBL(filtered)')
    drugs = rw.read_csv('./TheraNet/data/drugs(filtered)')
    pc_genes = rw.read_csv('./TheraNet/data/protein_coding_genes_ENSEMBL')
    PPI = rw.read_csv('./TheraNet/data/gene_gene_PPI700_ENSEMBL')
    PPI = PPI.fillna(0)
    PPI.loc[:,'max_ppi'] = PPI.loc[:,'max_ppi']/1000

    disease_genes = rw.read_csv('./TheraNet/data/'+disease_file)
    
    file = './TheraNet/output/drug_scores_TheraNet'
    candidate_drugs = rw.read_csv(file)
    
    network = matrix_to_network_TheraNet(drug_genes, PPI, disease_genes, pc_genes)
    
    scores = score_disease_genes(network,PPI, drug_genes, disease_genes, drugs, dgi, candidate_drugs)
    df = pd.DataFrame(scores)
    df.columns = disease_genes
    df = pd.concat([candidate_drugs, scores], ignore_index=True)
    rw.write_csv(df,'./TheraNet/output/drug_scores_TheraNet')
#..........................................................................................................
#
#
#TheraNetVSguney
#
#
#..........................................................................................................    
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3



def TheraNetVSguney():
    file = './TheraNet/output/drug_scores_guney'
    guney = rw.read_csv(file)

    file = './TheraNet/output/drug_scores_TheraNet'
    TheraNet = rw.read_csv(file)

    p_value = 0.05
    z_score = -1.96
    ind_zscore = [j for j, x in enumerate(TheraNet.loc[:,'z score']) if x  <= z_score]
    ind_pvalue = [j for j, x in enumerate(TheraNet.loc[:,'p-value']) if x  <= p_value]
    ind_TheraNet = intersection(ind_zscore, ind_pvalue)

    z_score = -0.15
    p_value = 1    
    ind_zscore = [j for j, x in enumerate(guney.loc[:,'z score']) if x  <= z_score]
    ind_pvalue = [j for j, x in enumerate(guney.loc[:,'p-value']) if x  <= p_value]
    ind_guney = intersection(ind_zscore, ind_pvalue)

    ind = intersection(ind_TheraNet, ind_guney)

    db = pd.DataFrame(index=range(1), columns=['TheraNet_method','guney_original method','mutually suggested drugs'])
    db.iloc[0,0] = len(ind_TheraNet)
    db.iloc[0,1] = len(ind_guney)
    db.iloc[0,2] = len(ind)
    print(f'{db}')

    TheraNet = TheraNet.loc[ind_TheraNet,:]
    TheraNet.to_csv('./TheraNet/output/suggested_drugs_TheraNet.csv', index=False)

    TheraNet = TheraNet.loc[ind_TheraNet,:]
    TheraNet.to_csv('./TheraNet/output/suggested_drugs_TheraNet.csv', index=False)
#..........................................................................................................
#
#
#pathways
#
#
#..........................................................................................................
def pathways():
    # Load datasets
    genes = pd.read_csv('./TheraNet/data/disease_target_genes.csv')  # Assumes a column "ENSEMBLE"
    pathways = pd.read_csv('./TheraNet/data/pathways_gene.csv')
    # Perform an efficient merge on the "ENSEMBLE" column
    merged_results = pathways.merge(genes, on="ENSEMBL", how="inner")
    # Save the result
    merged_results.to_csv('./TheraNet/output/affected_pathway_with_disease_target_genes.csv', index=False)    
#..........................................................................................................
#
#
#suggested_drugs
#
#
#..........................................................................................................  
def suggested_drugs(disease_file):
    #calling all the functions
    calculate_proximity_multiple(disease_file)
    check_for_drug_info()
    drug_scoring(disease_file)
    diease_target_gene_scoring(disease_file)
    TheraNetVSguney()
    pathways()
