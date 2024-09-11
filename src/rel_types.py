import time 
import networkx as nx
import numpy as np
import pandas as pd 
import sys 
import os 
import pickle
import argparse

import utils
from typing import List, Tuple, Dict, NamedTuple
import yaml 
import tqdm
from collections import namedtuple
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)




print("\tReading TPMS..")
tpms = pd.read_csv("../data/raw/expression/cri/iatlas-ici-hgnc_tpm.tsv",sep = "\t")
print("\tReading Counts...")
counts= pd.read_csv("../data/raw/expression/cri/iatlas-ici-hgnc_counts.tsv",sep = "\t")

counts['Run_ID'] = counts['Run_ID'].apply(lambda x: x.strip())
tpms['Run_ID'] = tpms['Run_ID'].apply(lambda x: x.strip())

tpm_genes = list(tpms.columns[1:])
count_genes = list(counts.columns[1:])

keep_genes = [x for x in count_genes if x in tpm_genes]

print(len(tpm_genes))
print(len(count_genes))
print(len(keep_genes))

STRING_aliases = pd.read_csv("../data/raw/networks/STRING/9606.protein.aliases.v12.0.txt",sep="\t")
STRING_aliases = STRING_aliases[STRING_aliases['alias'].isin(keep_genes)]
STRING_aliases = STRING_aliases[STRING_aliases['source'].isin(['BLAST_KEGG_NAME','BioMart_HUGO'])]
missing = [x for x in list(STRING_aliases['alias'].values) if x not in keep_genes]
print(missing)
print(len(missing))
