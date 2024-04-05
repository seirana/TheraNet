#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: s.hashemi

The function suggests drug for the set of disease target genes, 
        drug target the disease target genes, or they have PPI interaction with disease tare genes.

INPUT:
    ~/datat/disease_target_genes
    ~/data/drugs_per_gene
    ~/data/drugs_per_gene_PPI900
    
OUTPUT:
    ~/output/suggested_drugs
    ~/output/gene_drug_gene_matrix
"""

def suggested_drugs():

    import read_write as rw
    import pandas as pd
    
    file = '~/data/disease_target_genes'
    risk_genes = rw.read_csv(file)
    
    file = '~/data/drugs_per_gene'
    gene_drug = rw.read_csv(file)
       
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
    file = '~/data/drugs_per_gene_PPI900'
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
    
    file = '~/output/suggested_drugs'
    rw.write_csv(drugs,file)