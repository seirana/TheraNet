#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""
from IPython import get_ipython
get_ipython().run_line_magic('reset','-sf')

import os
os.system('clear')

import sys
sys.path.append('./TheraNet')

import TheraNet
disease_file = './TheraNet/data/disease_target_gene'
TheraNet.suggested_drugs(disease_file)
