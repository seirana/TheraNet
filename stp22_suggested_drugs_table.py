#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 12 2024

@author: s.hashemi
"""

from IPython import get_ipython
get_ipython().run_line_magic('reset','-sf')

import os
os.system('clear')

import sys
sys.path.append('/home/shashemi/Desktop/0-OneDrive/Code Python')

import read_write as rw
import pandas as pd

#direct targe genes
#file = '/home/shashemi/Desktop/0-OneDrive/Output/International PSC GWAS Study'
file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_regenie_genes'
risk_genes = rw.read_csv(file)

file = '/home/shashemi/Desktop/0-OneDrive/Data/DGIdb/drugs_per_gene'
gene_drug = rw.read_csv(file)

file = '/home/shashemi/Desktop/0-OneDrive/Data/STRING/gene_gene_PPI900'
ppi = rw.read_csv(file)

#identify drugs
drugs = pd.DataFrame(index=range(len(risk_genes)),columns=['drugs','drugs_PPI900','drug_target_genes'])
for i in range(len(risk_genes)):
    indices = [j for j, x in enumerate(gene_drug.loc[:,'ENSEMBL']) if x == risk_genes.loc[i,'ENSEMBL']]
    if len(indices) > 0:
        tmp = ''
        for j in range(len(indices)):
            tmp = tmp + '#' +  gene_drug.loc[indices[j],'drug_claim_name']
        drugs.loc[i,'drugs'] =  tmp

#identify drugs with PPI900
file = '/home/shashemi/Desktop/0-OneDrive/Data/DGIdb/drugs_per_gene_PPI900'
gene_drug = rw.read_csv(file)

for i in range(len(risk_genes)):
    indices = [j for j, x in enumerate(gene_drug.loc[:,'genes_PPI1']) if risk_genes.loc[i,'ENSEMBL'] in x]
    if len(indices) > 0:
        tmp = ''
        for j in range(len(indices)):
            tmp = tmp + '#' +  gene_drug.loc[indices[j],'drug_claim_name']
        drugs.loc[i,'drugs_PPI900'] =  tmp

        tmp = ''
        for j in range(len(indices)):
            tmp = tmp + '#' +  gene_drug.loc[indices[j],'Gene stable ID']
        drugs.loc[i,'drug_target_genes'] =  tmp

file = '/home/shashemi/Desktop/0-OneDrive/Output/Vis/PSC_suggested_drugs(drug_network)'
rw.write_csv(drugs,file)