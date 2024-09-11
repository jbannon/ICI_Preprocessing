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

""" 
Code to take in the raw data from the CRI iAtlas and:

1)  output the .csv files for response and TPM measurements by tissue
2)  compute and write the DE genesets
"""







def main():

	drugs = ['Atezo','Pembro','Nivo','Ipi','Ipi + Pembro', 'Ipi + Nivo']
	# drugs = ['Ipi']
	# _valid_drugs = ['Atezo','Nivo','Pembro','Ipi','Ipi + Pembro']
	VALID_STRING_SOURCES = ['BLAST_KEGG_NAME','BioMart_HUGO']
	# TIDE	TIDE_IFNG	TIDE_MSI	TIDE_CD274	TIDE_CD8	TIDE_CTL_Flag	TIDE_Dysfunction	TIDE_Exclusion	TIDE_MDSC	TIDE_CAF	TIDE_TAM_M2	TIDE_CTL
	n_iters = 200
	print("\n\tReading TPM...")
	tpms = pd.read_csv("../data/raw/expression/cri/iatlas-ici-hgnc_tpm.tsv",sep = "\t")
	print("\tReading Counts...")
	counts= pd.read_csv("../data/raw/expression/cri/iatlas-ici-hgnc_counts.tsv",sep = "\t")
	print("\tReading Features..")
	feats= pd.read_csv("../data/raw/expression/cri/iatlas-ici-features.tsv",sep = "\t")
	
	kept_immune_features = ['IMPRES','Miracle','TIDE','TIDE_IFNG',
		'TIDE_MSI','TIDE_CD274','TIDE_CD8','TIDE_CTL_Flag',
		'TIDE_Dysfunction', 'TIDE_Exclusion','TIDE_MDSC',
		'TIDE_CAF','TIDE_TAM_M2','TIDE_CTL']
	
	feats = feats[['Run_ID']+kept_immune_features]
	

	print("\tReading Response Data..")
	response_data = pd.read_csv("../data/raw/expression/cri/iatlas-ici-sample_info.tsv", sep = "\t")
	response_data = response_data.dropna(subset=['Response'])
	response_data['Binned_Response'] = response_data['Response'].apply(lambda x: "R" if x in ['Complete Response','Partial Response'] else "NR")
	response_data['Response'] = pd.Categorical(response_data['Binned_Response'])
	response_data['Binned_Response'] = pd.Categorical(response_data['Binned_Response'])
	response_data['Response'] = response_data['Response'].cat.codes
	response_data = response_data[response_data['Sample_Treated'] == False]
	response_data = response_data[response_data['ICI_Tx'] == True]
	

	
	response_data = response_data[['Run_ID','TCGA_Tissue','ICI_Rx','Response','Binned_Response']]

	response_data['Run_ID'] = response_data['Run_ID'].apply(lambda x: x.strip())
	feats['Run_ID'] = feats['Run_ID'].apply(lambda x: x.strip())
	counts['Run_ID'] = counts['Run_ID'].apply(lambda x: x.strip())
	tpms['Run_ID'] = tpms['Run_ID'].apply(lambda x: x.strip())
	

	STRING_aliases = pd.read_csv("../data/raw/networks/STRING/9606.protein.aliases.v12.0.txt",sep="\t")
	STRING_aliases = STRING_aliases[STRING_aliases['source'].isin(VALID_STRING_SOURCES)]
	STRING_links = pd.read_csv("../data/raw/networks/STRING/9606.protein.links.v12.0.txt",sep = ' ')
	
	

	tpm_genes = set(tpms.columns[1:])
	count_genes = set(counts.columns[1:])
	measured_genes = [x for x in tpm_genes if x in count_genes]
	

	tpms_full = tpms.copy(deep=True)
	counts_full = counts.copy(deep = True)
	tpms_full = tpms_full[['Run_ID']+measured_genes]
	counts_full = counts_full[['Run_ID']+measured_genes]
	good_pids, keep_genes, pid_gene_map = utils.find_unambiguous_nodes(STRING_aliases,measured_genes)
	
	
	tpms = tpms[['Run_ID']+list(keep_genes)]
	counts = counts[['Run_ID']  + list(keep_genes)]
	# print(tpms.columns)
	
	
	print("\tMaking Networks...")
	STRING_links = STRING_links[(STRING_links['protein1'].isin(good_pids)) & (STRING_links['protein2'].isin(good_pids))]
	

	sparse, sparseW, dense, denseW, easy, easyW, tight, tightW=\
		utils.construct_networks(STRING_links,pid_gene_map, fname = "../data/processed/networks/cri/stats_.txt")


	
	netList = [sparse,sparseW,dense, denseW,easy, easyW,tight,tightW]
	for network, name,stem in tqdm.tqdm(zip(netList,
		['unweighted','weighted','unweighted','weighted','unweighted','weighted','unweighted','weighted'], 
		["sparse","sparse","dense","dense", "easy", "easy","tight", "tight"]),total = len(netList)):
		
		os.makedirs("../data/processed/networks/cri/{s}".format(s=stem),exist_ok = True)
		with open("../data/processed/networks/cri/{s}/{n}.pickle".format(s = stem,n =name),'wb') as ostream:
			pickle.dump(network,ostream)


	print("\n\tProcessing Expression...")
	
	for drug in drugs: 
		treated_samples = response_data[response_data['ICI_Rx']==drug]
		tissues = pd.unique(treated_samples['TCGA_Tissue'])
		if drug == "Ipi + Pembro":
			drug = "Ipi+Pembro"
		elif drug == "Ipi + Nivo":
			drug = "Ipi+Pembro"
		
		for tissue in tissues:
			
			print("\n\tWorking on {d}, {t}".format(d=drug,t=tissue))


			os.makedirs("../data/processed/expression/cri/{d}/{t}/".format(d=drug,t=tissue),exist_ok = True)

			responses = treated_samples[treated_samples['TCGA_Tissue']==tissue][['Run_ID','Response']]
			print("IDs:")
			print(responses['Run_ID'])
			continue
			expression_measurements = tpms[tpms['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			
			
			full_expression_measurements = tpms_full[tpms_full['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			immune_features = feats[feats['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			counts_data = counts[counts['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			full_count_data = counts_full[counts_full['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			
			
			responses.set_index('Run_ID',inplace = True)
			sample_order = list(responses.index)
			
			responses.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			full_expression_measurements.set_index('Run_ID',inplace = True)
			full_expression_measurements = full_expression_measurements.loc[sample_order,:]
			full_expression_measurements.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			full_count_data.set_index('Run_ID',inplace = True)
			full_count_data = full_count_data.loc[sample_order,:]
			full_count_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			expression_measurements.set_index('Run_ID',inplace = True)
			expression_measurements = expression_measurements.loc[sample_order,:]
			expression_measurements.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			counts_data.set_index('Run_ID',inplace = True)
			counts_data = counts_data.loc[sample_order,:]
			counts_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])


			immune_features.set_index('Run_ID',inplace = True)
			immune_features = immune_features.loc[sample_order,:]
			immune_features.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			print("\n")
			print("Full Expression:\t{s}".format(s=full_expression_measurements.shape))
			print("Full Count:\t{s}".format(s=full_count_data.shape))
			print("Expression TPM:\t{s}".format(s=expression_measurements.shape))
			print("Expression Counts:\t{s}".format(s=counts_data.shape))
			print("Immune Features:\t{s}".format(s=immune_features.shape))
			print("\n")
		
			full_expression_measurements.to_csv("../data/processed/expression/cri/{d}/{t}/expression_full.csv".format(d=drug,t=tissue),index = False)
			full_count_data.to_csv("../data/processed/expression/cri/{d}/{t}/counts_full.csv".format(d=drug,t=tissue),index = False)
			
			expression_measurements.to_csv("../data/processed/expression/cri/{d}/{t}/expression.csv".format(d=drug,t=tissue),index = False)
			responses.to_csv("../data/processed/expression/cri/{d}/{t}/response.csv".format(d=drug,t=tissue),index = False)
			immune_features.to_csv("../data/processed/expression/cri/{d}/{t}/features.csv".format(d=drug,t=tissue),index = False)
			counts_data.to_csv("../data/processed/expression/cri/{d}/{t}/counts.csv".format(d=drug,t=tissue),index = False)
			

			print("\tcomputing differentially expressed genes...")
			if tissue == 'GBM':
				continue
			col_data = treated_samples[treated_samples['TCGA_Tissue']==tissue][['Run_ID','Binned_Response']]
			
			count_matrix = counts[counts['Run_ID'].isin(list(col_data['Run_ID']))]
			
			col_data.reset_index(inplace = True, drop = True)
			count_matrix.reset_index(inplace = True, drop = True)
			col_data.columns = ['Run_ID','Binned_Response']
			col_data.set_index('Run_ID',inplace = True)
			sample_order = list(col_data.index)
			col_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			col_data.to_csv("./col_data.csv",index = False)
			
			count_matrix.set_index('Run_ID',inplace = True)
			count_matrix = count_matrix.loc[sample_order,:]
			count_matrix.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			count_matrix = count_matrix.transpose()
			count_matrix.reset_index(inplace = True)
			count_matrix.columns = count_matrix.iloc[0]
			count_matrix.drop(index=0,axis=0,inplace=True)
			count_matrix.rename(columns = {'Run_ID':'Gene_ID'},inplace = True)
			count_matrix.to_csv("./count_data.csv",index = False)
			
			
			
			cmd = "Rscript compute_DE_genes.R {drug} {tissue} {n_iters} ./count_data.csv ./col_data.csv".format(drug=drug, 
				tissue = tissue, 
				n_iters = n_iters)
			os.system(cmd)

			os.remove("./col_data.csv")
			os.remove("./count_data.csv")
		

		if drug in ['Atezo','Pembro','Nivo']:
			print("\n\tWorking on {d}, PANCAN".format(d=drug))
			os.makedirs("../data/processed/expression/cri/{d}/PANCAN/".format(d=drug),exist_ok = True)
			responses = treated_samples[['Run_ID','Response']]
			

			expression_measurements = tpms[tpms['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			immune_features = feats[feats['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			counts_data = counts[counts['Run_ID'].isin(pd.unique(responses['Run_ID']))]

			responses.set_index('Run_ID',inplace = True)
			sample_order = list(responses.index)
			
			responses.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			full_expression_measurements = tpms_full[tpms_full['Run_ID'].isin(pd.unique(responses['Run_ID']))]
			full_count_data = counts_full[counts_full['Run_ID'].isin(pd.unique(responses['Run_ID']))]

			full_expression_measurements.set_index('Run_ID',inplace = True)
			full_expression_measurements = full_expression_measurements.loc[sample_order,:]
			full_expression_measurements.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			full_count_data.set_index('Run_ID',inplace = True)
			full_count_data = full_count_data.loc[sample_order,:]
			full_count_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			expression_measurements.set_index('Run_ID',inplace = True)
			expression_measurements = expression_measurements.loc[sample_order,:]
			expression_measurements.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			
			counts_data.set_index('Run_ID',inplace = True)
			counts_data = counts_data.loc[sample_order,:]
			counts_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])


			immune_features.set_index('Run_ID',inplace = True)
			immune_features = immune_features.loc[sample_order,:]
			immune_features.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			print("\n")
			print("Full Expression:\t{s}".format(s=full_expression_measurements.shape))
			print("Full Count:\t{s}".format(s=full_count_data.shape))
			print("Expression TPM:\t{s}".format(s=expression_measurements.shape))
			print("Expression Counts:\t{s}".format(s=counts_data.shape))
			print("Immune Features:\t{s}".format(s=immune_features.shape))
			print("\n")

			full_expression_measurements.to_csv("../data/processed/expression/cri/{d}/PANCAN/expression_full.csv".format(d=drug,t=tissue),index = False)
			full_count_data.to_csv("../data/processed/expression/cri/{d}/PANCAN/counts_full.csv".format(d=drug,t=tissue),index = False)
			expression_measurements.to_csv("../data/processed/expression/cri/{d}/PANCAN/expression.csv".format(d=drug),index = False)
			responses.to_csv("../data/processed/expression/cri/{d}/PANCAN/response.csv".format(d=drug),index = False)
			immune_features.to_csv("../data/processed/expression/cri/{d}/PANCAN/features.csv".format(d=drug),index = False)
			counts_data.to_csv("../data/processed/expression/cri/{d}/PANCAN/counts.csv".format(d=drug),index = False)
			

			print("\tcomputing differentially expressed genes...")
			col_data = treated_samples[['Run_ID','Binned_Response']]
			
			count_matrix = counts[counts['Run_ID'].isin(list(col_data['Run_ID']))]
			
			
			col_data.reset_index(inplace = True, drop = True)
			count_matrix.reset_index(inplace = True, drop = True)
			col_data.columns = ['Run_ID','Binned_Response']
			col_data.set_index('Run_ID',inplace = True)
			sample_order = list(col_data.index)
			col_data.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			col_data.to_csv("./col_data.csv",index = False)
			
			count_matrix.set_index('Run_ID',inplace = True)
			count_matrix = count_matrix.loc[sample_order,:]
			count_matrix.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			count_matrix = count_matrix.transpose()
			count_matrix.reset_index(inplace = True)
			count_matrix.columns = count_matrix.iloc[0]
			count_matrix.drop(index=0,axis=0,inplace=True)
			count_matrix.rename(columns = {'Run_ID':'Gene_ID'},inplace = True)
			count_matrix.to_csv("./count_data.csv",index = False)
			
			

			cmd = "Rscript compute_DE_genes.R {drug} {tissue} {n_iters} ./count_data.csv ./col_data.csv".format(drug=drug, 
				tissue = "PANCAN", 
				n_iters = n_iters)
			os.system(cmd)
			os.remove("./col_data.csv")
			os.remove("./count_data.csv")





		
	
	








if __name__ == '__main__':
	main()