import numpy as np
from collections import defaultdict
import scipy.stats as stat
import pandas as pd
import time, os, random, sys, argparse
import networkx as nx
from networkx.algorithms.link_analysis import pagerank
import tqdm as tqdm
import yaml 
from typing import List, Tuple, Dict, NamedTuple
from collections import namedtuple, defaultdict

def fetch_gdsc_targets(
	chebi_dict: Dict[str,str],
	drug_info: pd.DataFrame
	) -> Dict[str,list[str]]:
	tgts = {}
	for drug in tqdm.tqdm(chebi_dict.keys()):	
		cid = chebi_dict[drug]
		temp = drug_info[drug_info['PARTICIPANT_A']==cid]
		tgts[drug] = list(pd.unique(temp['PARTICIPANT_B']))

	return tgts


def build_geneset(
	fpath:str,
	lincs:bool = False
	) -> None:
	seperator = "\t" if lincs else ","
	symcol = 'Symbol' if lincs else 'Gene Symbol'
	df = pd.read_csv(fpath,sep=seperator)
	
	genes = list(pd.unique(df[symcol]))
	return genes
	
def get_hallmark_genes(
    file_path:str
    )->List[str]:
    
    gene_set = []
    with open(file_path,"r") as f:
        lines = f.readlines()

    for line in lines:
        temp = line.split("\t")
        temp = [x.rstrip() for x in temp ]
        temp = temp[2:]
        gene_set = gene_set + temp

    gene_set = list(set(gene_set))
    return gene_set

def make_chebi_dict(
	fname:str
	)-> Dict[str,str]:
	with open(fname,'r') as istream:
		lines = istream.readlines()
	lines = [line.rstrip() for line in lines]
	lines = [line.split("\t") for line in lines]
	drug_to_chebi = {line[0]:line[1] for line in lines[1:]}	
	return drug_to_chebi


def fetch_targets(
	drug_info: pd.DataFrame,
	chebi_dict: Dict[str,str],
	ICI_targets:Dict[str,List[str]] = {'PD-L1':['CD274'],'PD1':['PDCD1'],'CTLA4':['CTLA4'], 'CTLA4+PD1': ['CTLA4','PDCD1']}
	) -> Dict[str,List[str]]:
	tgts = fetch_gdsc_targets(chebi_dict,drug_info)
	tgts.update(ICI_targets)
	return tgts
	

def main(config:Dict):

	hmg = get_hallmark_genes(config['hallmark_path'])
	with open("../data/processed/genesets/mdsig_hallmarks.txt", "w") as f:
		f.writelines([g+"\n" for g in hmg])
	lincs_gs = build_geneset(config['lincs_path'],True)
	with open("../data/processed/genesets/LINCS.txt", "w") as f:
		f.writelines([g+"\n" for g in lincs_gs])
	cosmic_gs = build_geneset(config['cosmic_path'])
	with open("../data/processed/genesets/COSMIC.txt", "w") as f:
		f.writelines([g+"\n" for g in cosmic_gs])
	
	auslander_map = {'PD-L1':'CD274','PD1':'PDCD1','CTLA4':'CTLA4','PDL-1':'CD274','PDL1':'CD274',
					'GAL3':'LGAL3','OX-40L':'TNFSF4','OX40L':'TNFSF4', 'TIM-3':'HAVCR2','PD-1LG2':'PDCD1LG2','CD70LG':'CD70',
					'VISTA':'VSIR','CD266':'TNFRSF12A', 'CD40L':'CD40LG','DR3':'TNFRSF25','ICOSL':'ICOSLG','NAIL':'CD244',
					'PD-1':'PDCD1','PVRL2':'NECTIN2', 'SLAM':'SLAMF1','TIM2':'HAVCR2','HVEM':'TNFRSF14','CD137L':'TNFSF9',
					'LGAL3':'LGALS3'}
					
	print("loading data + mapping targets")
	aliases = pd.read_csv(config['alias_path'], sep=config['alias_sep']) #corresponds to netprop annotation
	links = pd.read_csv(config['link_path'], sep=config['link_sep']) 
	string_cutoff = config['string_cutoff']
	drug_chebi = make_chebi_dict(config['drug_chebi_path'])
	drugbank_info = pd.read_csv(config['drugbank_path'],sep=config['drugbank_sep'])
	drugbank_info = drugbank_info[drugbank_info['INTERACTION_TYPE']=='chemical-affects']
	targets = fetch_targets(drugbank_info,drug_chebi)

	print("harmonizing auslander genes")

	with open("../data/raw/genesets/auslander_45.txt","r") as f:
		lines = f.readlines()

	auslander_genes = [line.rstrip() for line in lines]

	cri_tpm = pd.read_csv(config['cri_path'],sep=config['cri_sep'],nrows = config['cri_rows'])
	tpm_cols  = list(set(cri_tpm.columns[1:]))
	
	ag_ = []
	for gene in auslander_genes:
		if gene in auslander_map.keys():
			gene = auslander_map[gene]
		if gene not in tpm_cols:
			print("{g} not in".format(g=gene))
		else:
			ag_.append(gene)
	
	with open("../data/processed/genesets/auslander.txt", "w") as f:
		f.writelines([g+"\n" for g in ag_])

	
	print("\n")
	print("building network")

	
	tmp_G = nx.Graph()

	nodes1 = links.values[:,0]
	nodes2 = links.values[:,1]
	scores = links.values[:,2]
	
	for n1, n2, score in tqdm.tqdm(zip(nodes1, nodes2, scores),total=len(nodes1)):
		if score >= string_cutoff:
			tmp_G.add_edge(n1, n2)
	
	LCC_genes = max(nx.connected_components(tmp_G), key=len)

	G = tmp_G.subgraph(LCC_genes) 
	network_nodes = G.nodes()
	network_edges = G.edges()

	print('network nodes: %s'%len(network_nodes)) # 16584
	print('network edges: %s'%len(network_edges))
	print('\n')

	anno_dic = defaultdict(list) 
	anno = aliases[aliases['source'].isin(config['valid_sources'])] # my own selection
	for idx, row in anno.iterrows():
		anno_dic[row['#string_protein_id']].append(row['alias'])

	path = "../data/processed/genesets/"
	os.makedirs(path,exist_ok=True)

	for target in targets.keys():
		print('\ttesting %s, %s'%(target, time.ctime()))
		print('\tcorresponding targets: {g}'.format(g=targets[target]))

		output = defaultdict(list)
		output_col = ['gene_id', 'string_protein_id', 'propagate_score']

		# network propagation
		biomarker_genes = targets[target]
		
		pIDs = aliases.loc[aliases['alias'].isin(biomarker_genes),:]['#string_protein_id'].tolist() #only remain biomarker_gene

		propagate_input = {}
		for node in network_nodes:
			propagate_input[node] = 0
			if node in pIDs:
				propagate_input[node] = 1
		
		propagate_scores = pagerank(G, personalization=propagate_input, max_iter=100, tol=1e-06) ## NETWORK PROPAGATION
		
		for prot in propagate_scores.keys():
			if prot in list(anno_dic.keys()):
				for g in anno_dic[prot]:
					output['geneID'].append(g)
					output['string_protein_id'].append(prot)
					output['propagate_scores'].append(propagate_scores[prot])
		
		output = pd.DataFrame(data=output)
		output = output.sort_values(by=['propagate_scores'], ascending=False)
		output.to_csv("{p}{t}.txt".format(p=path,t=target),sep='\t',index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-config",help="The config file for these experiments")
	args = parser.parse_args()

	with open(args.config) as file:
		config = yaml.safe_load(file)

	main(config)
		