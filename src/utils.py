import networkx as nx 
import tqdm 
import sys
from typing import List, Tuple, Dict, NamedTuple
import numpy as np 
import pandas as pd 

def find_unambiguous_nodes(
	aliases:pd.DataFrame,
	measured_genes:List[str]
	) -> Tuple[List[str],List[str]]:

	common_genes =  list(set(measured_genes).intersection(set(aliases['alias'])))	
	ambiguous_genes = []
	good_genes = []
	good_string_ids = []
	empty_genes = []
	protein_gene_map = {}
	print("\n reconciling measured genes and STRING aliases\n")

	for gene in tqdm.tqdm(common_genes):
		temp = aliases[aliases['alias'] == gene]

		if len(pd.unique(temp['#string_protein_id']))>1:
			ambiguous_genes.append(gene)
		elif len(pd.unique(temp['#string_protein_id']))==1:
			good_genes.append(gene)
			# print(pd.unique(temp['#string_protein_id']))
			good_string_ids.append( list(pd.unique(temp['#string_protein_id']))[0])
			protein_gene_map[list(pd.unique(temp['#string_protein_id']))[0]] = gene
		elif temp.shape[0]==0:
			empty_genes.append(gene)
			print('empty')

	stat_string = "\n\t Started with {cg} genes\n\t Removed {ag} ambiguous genes\n\t {mg} missing genes\n\t keeping {gg}".format(
		cg = len(common_genes),
		ag = len(ambiguous_genes),
		mg = len(empty_genes),
		gg = len(good_genes))
	
	print(stat_string)
	
	

	return good_string_ids, good_genes, protein_gene_map
	
def construct_networks(
	network:pd.DataFrame,
	protein_gene_map:Dict[str,str],
	sparse_cutoff:int = 850,
	tight_cutoff:int = 950,
	fine_cutoff:int = 990,
	dense_cutoff:int = 600,
	easy_cutoff:int = 500,
	fname:str = None
	) -> Tuple[nx.Graph, nx.Graph, nx.Graph, nx.Graph]:
	nodes1 = network['protein1'].values
	nodes2 = network['protein2'].values
	scores = network['combined_score'].values
	
	sparse = nx.Graph()
	sparse_weighted = nx.Graph()

	fine = nx.Graph()
	fine_weighted = nx.Graph()

	dense = nx.Graph()
	dense_weighted = nx.Graph()

	easy = nx.Graph()
	easy_weighted = nx.Graph()

	tight = nx.Graph()
	tight_weighted = nx.Graph()
	
	print('constructing networks')
	for n1, n2, score in tqdm.tqdm(zip(nodes1, nodes2, scores),total = network.shape[0]):
		# print(row)
		# n1, n2, score = row['protein1'], row['protein2'], row['combined_score']
		node1 = protein_gene_map[n1]
		node2 = protein_gene_map[n2]
		if score >= fine_cutoff:
			fine.add_edge(node1, node2)
			fine_weighted.add_edge(node1,node2,weight=score)
		if score >= tight_cutoff:
			tight.add_edge(node1, node2)
			tight_weighted.add_edge(node1,node2,weight=score)
		if score >= sparse_cutoff:
			sparse.add_edge(node1, node2)
			sparse_weighted.add_edge(node1,node2,weight=score)
		if score >= dense_cutoff:
			dense.add_edge(node1, node2)
			dense_weighted.add_edge(node1,node2, weight=score)
		if score >= easy_cutoff:
			easy.add_edge(node1, node2)
			easy_weighted.add_edge(node1,node2, weight=score)
	

	LCC_genes = max(nx.connected_components(sparse), key=len)
	network_string = "".join(["Sparse Network Has:\n", "\t {n} nodes\n".format(n=len(sparse.nodes())),
		 "\t {e} edges\n".format(e=len(sparse.edges())),
		  "\t {n}-node LCC\n".format(n=len(LCC_genes))])


	LCC_genes = max(nx.connected_components(fine), key=len)
	network_string = network_string +"".join(["Fine Network Has:\n",
		"\t {n} nodes\n".format(n=len(fine.nodes())),
		"\t {e} edges\n".format(e=len(fine.edges())),
		"\t {n}-node LCC\n".format(n=len(LCC_genes))])

	LCC_genes = max(nx.connected_components(dense), key=len)
	network_string = network_string +"".join(["Dense Network Has:\n",
		"\t {n} nodes\n".format(n=len(dense.nodes())),
		"\t {e} edges\n".format(e=len(dense.edges())),
		"\t {n}-node LCC\n".format(n=len(LCC_genes))])

	LCC_genes = max(nx.connected_components(easy), key=len)
	network_string = network_string +"".join(["Easy Network Has:\n",
		"\t {n} nodes\n".format(n=len(easy.nodes())),
		"\t {e} edges\n".format(e=len(easy.edges())),
		"\t {n}-node LCC\n".format(n=len(LCC_genes))])

	LCC_genes = max(nx.connected_components(tight), key=len)
	network_string = network_string +"".join(["Tight Network Has:\n", "\t {n} nodes\n".format(n=len(tight.nodes())),
		 "\t {e} edges\n".format(e=len(tight.edges())),
		  "\t {n}-node LCC\n".format(n=len(LCC_genes))])

	with open(fname,"w") as ostream:
		ostream.write(network_string)
	
	return sparse, sparse_weighted, dense, dense_weighted, easy, easy_weighted, tight, tight_weighted
