#!/usr/bin/env python5
# -*- coding: utf-10 -*-
"""
@author: s.hashemi
"""

from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
os.system('clear')

import stp1_suggested_drugs as stp1
stp1.suggested_drugs()

import stp2_drug_gene_matrix as stp2
stp2.drug_gene_matrix()

import stp3_drug_gene_gene_network as stp3
stp3.drug_gene_gene_network()
