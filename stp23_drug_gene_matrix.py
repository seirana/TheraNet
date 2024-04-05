#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:48:31 2024

@author: shashemi
"""

from IPython import get_ipython
get_ipython().run_line_magic('reset','-sf')

import os
os.system('clear')

import sys
sys.path.append('/home/shashemi/Desktop/0-OneDrive/Code Python')

import read_write as rw
import pandas as pd
import numpy as np

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def convert(string):
    string = string.replace('##','')
    if isinstance(string, str) == False:
        string = string.to_string()
    lst = list(string.split('#'))
    return lst

file = '/home/shashemi/Desktop/0-OneDrive/Data/BioMart/gene'
gene_id = rw.read_csv(file)

file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_suggested_drugs(drug_network)'
sug_drug = rw.read_csv(file)

for i in range(len(sug_drug)):
    d = sug_drug.loc[i,'drugs']
    if isfloat(d) == True:
        sug_drug.loc[i,'drugs'] = 'nan'
        
    d = sug_drug.loc[i,'drugs_PPI900']
    if isfloat(d) == True:
        sug_drug.loc[i,'drugs_PPI900'] = 'nan'

    d = sug_drug.loc[i,'drug_target_genes']
    if isfloat(d) == True:
        sug_drug.loc[i,'drug_target_genes'] = 'nan'
        
for i in range(len(sug_drug)):
    a = sug_drug.loc[i,'drugs']
    sug_drug.loc[i,'drugs'] = a.upper()
    a = sug_drug.loc[i,'drugs_PPI900']
    sug_drug.loc[i,'drugs_PPI900'] = a.upper()
    a = sug_drug.loc[i,'drug_target_genes']
    sug_drug.loc[i,'drug_target_genes'] = a.upper()
        
#direct targe genes
#file = '/home/shashemi/Desktop/0-OneDrive/Output/International PSC GWAS Study'
file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_regenie_genes'
risk_gene = rw.read_csv(file)

drug = []
for i in range(len(sug_drug)):
    d = sug_drug.loc[i,'drugs']
    if isfloat(d) == False:
        l = convert(d)
        drug = drug + l       

    d = sug_drug.loc[i,'drugs_PPI900']
    if isfloat(d) == False:
        l = convert(d)
        drug = drug + l

drug = pd.DataFrame(drug)
drug = drug.drop_duplicates()
drug = drug.drop(drug.index[0])
drug.index = range(len(drug))
drug.columns = ['drugs']

file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_nodes_drug_list(drug_network)'
rw.write_csv(drug,file)

gene = []
for i in range(len(risk_gene)):
    g = sug_drug.loc[i,'drugs_PPI900']
    if isfloat(g) == False:
        gl = convert(sug_drug.loc[i,'drug_target_genes'])
        gene = gene + gl

gene = pd.DataFrame(gene)
gene = gene.drop_duplicates()
gene.index = range(len(gene))

gene = pd.concat([risk_gene.loc[:,'ENSEMBL'],gene],axis=0)
gene = gene.drop_duplicates()
gene.columns = ['genes']
gene = gene.sort_values(by='genes')
gene.index = range(len(gene))
l = len(gene)
gene = gene.loc[1:l,'genes']
gene.index = range(len(gene))
gene = pd.DataFrame(gene,columns=['genes'])

gene_type = pd.DataFrame(index=range(len(gene)),columns=['gene_id','gene_type'])
for i in range(len(gene)):
    
    ind = [j for j, x in enumerate(gene_id.loc[:,'Gene stable ID']) if x == gene.loc[i,'genes']]
    if len(ind)>0:
        gene_type.loc[i,'gene_id'] = gene_id.loc[ind[0],'Gene name']
    else:
        ind = [j for j, x in enumerate(risk_gene.loc[:,'ENSEMBL']) if x == gene.loc[i,'genes']]
        gene_type.loc[i,'gene_id'] = risk_gene.loc[ind[0],'ENSEMBL']
        
    ind = [j for j, x in enumerate(risk_gene.loc[:,'ENSEMBL']) if x == gene.loc[i,'genes']]
    if len(ind) > 0:
        g = sug_drug.loc[ind,'drugs']
        if isfloat(g) == False:
            gene_type.loc[i,'gene_type'] = 'drug_and_disease_target'
        else:
            g = sug_drug.loc[ind,'drugs_PPI900']
            if isfloat(g) == False:
                gene_type.loc[i,'gene_type'] = 'indirect_drug_target'
            else:
                gene_type.loc[i,'gene_type'] = 'non_drug_target'
    else:
        gene_type.loc[i,'gene_type'] = 'drug_target'

db = pd.concat([gene,gene_type],axis =1)
file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_nodes_gene_list_type(drug_network)'
rw.write_csv(db,file)

mat = np.zeros((len(gene),len(drug)+len(gene)), dtype=int)
for i in range(len(gene)): #gene-drug-edges
    ind = [j for j, x in enumerate(risk_gene.loc[:,'ENSEMBL']) if x == gene.loc[i,'genes']]
    if len(ind) > 0: #risk_gene
        d = sug_drug.loc[ind,'drugs'] 
        if isfloat(d) == False: #direct drug
            d = convert(d)
            d = pd.DataFrame(d)
            d = d.drop_duplicates()
            d.index = range(len(d))

            for k in range(len(d)):
                ind1 = [j for j, x in enumerate(drug.loc[:,'drugs']) if x == d.loc[k,0]]
                if len(ind1)>0:
                    mat[i,ind1] = 1
    else: #indirect drug
        ind = [j for j, x in enumerate(sug_drug.loc[:,'drug_target_genes']) if gene.loc[i,'genes'] in x]
        if len(ind)>1:
            dt = str()
            for k in range(len(ind)):
                a = sug_drug.loc[ind[k],'drugs_PPI900']
                dt = dt+'#'+a
                d = dt
        else:
            d = sug_drug.loc[ind[0],'drugs_PPI900']
                
        if isfloat(d) == False:
            d = convert(d)
            d = pd.DataFrame(d)
            d = d.drop_duplicates()
            d.index = range(len(d))

            for k in range(len(d)):
               ind1 = [j for j, x in enumerate(drug.loc[:,'drugs']) if x == d.loc[k,0]]
               if len(ind1)>0:
                   mat[i,ind1] = 1
        

file = '/home/shashemi/Desktop/0-OneDrive/Data/STRING/gene_gene_PPI900'
ppi = rw.read_csv(file)

for i in range(len(gene)): #gene-gene-edges 
    ind1 = [j for j, x in enumerate(ppi.loc[:,'gene1']) if x == gene.loc[i,'genes']]
    if len(ind1) > 0:
        for k in range(len(ind1)):
            ind2 = [j for j, x in enumerate(gene.loc[:,'genes']) if x == ppi.loc[ind1[k],'gene2']]
            if len(ind2)>0:
                mat[i,len(drug)+ind2[0]] = 1
                mat[ind2[0],len(drug)+i] = 1

mat = pd.DataFrame(mat)
file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_gene_drug_gene_matrix(drug_network)'
rw.write_csv(mat,file)