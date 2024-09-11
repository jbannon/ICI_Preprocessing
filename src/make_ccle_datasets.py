import time 
import networkx as nx
import numpy as np
import pandas as pd 
import sys 
import os 
import pickle
import argparse
import matplotlib.pyplot as plt 

import utils
from typing import List, Tuple, Dict, NamedTuple
import yaml 
import tqdm
from collections import namedtuple
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

def squeezeString(
	s:str,
	delims = [" ","-",".","/", ";", ":", "(",")", ","]
	)->str:
	for delim in delims:
		s = "".join(s.split(delim))
	s = s.upper()
	return s




def fetch_ccl_annotations(
	annFile:str,
	cclMetaFile:str
	):
	

	

	cdf1 = pd.read_csv(annFile,sep = "\t") ## TCGA source
	cdf1 = cdf1[['CCLE_ID','Name','tcga_code']]
	print(cdf1['tcga_code'].value_counts())
	print(np.sum(cdf1['tcga_code'].value_counts()))
	sys.exit(1)
	cdf1.dropna(inplace = True)
	cdf1['Name'] = cdf1["Name"].map(squeezeString)
	
	cdf2 = pd.read_csv(cclMetaFile,sep="\t")
	cdf2 = cdf2[['master_ccl_id','ccl_name']]

	cdf2.columns = ['master_ccl_id','Name']
	
	cdf2["Name"] = cdf2['Name'].map(squeezeString)

	cdf = cdf1.merge(cdf2,on="Name")

	return cdf


def fetch_compound_info(
	compoundFile:str, 
	assayFile:str,
	responseFile:str,
	drugs:List[str]
	):


	compound = pd.read_csv(compoundFile,sep = "\t")
	compound = compound[['master_cpd_id','cpd_name','gene_symbol_of_protein_target']]
	compound = compound[compound['cpd_name'].isin(drugs)]
	assay = pd.read_csv(assayFile,sep="\t")
	assay = assay[['experiment_id','master_ccl_id']]
	respInfo = pd.read_csv(responseFile,sep = "\t")
	respInfo = respInfo[['experiment_id','master_cpd_id','area_under_curve']]
	return compound, assay, respInfo

def main():

	"""
	Crizotinib: {source:'LUAD','BRCA'}
	Sunitinib: {source:'LUAD','SKCM'}
	Sorafenib: {source: 'LUAD','BRCA','SKCM'} target: LIHC, LUSC
	Erlotinib: {source: 'LUAD','BRCA',} target: 'OV','BLCA', 'LIHC','LUSC'

	"""
	tissueDict = {'LUAD':['LUAD'],
		'BRCA':['BRCA'],
		'SKCM':['SKCM'],
		'BROAD':['LUAD','COAD/READ','SKCM','PAAD','OV','PAAD','KIRC','ESCA','LUSC','HNSC','BRCA','STAD','LUSC','LIHC']}
	
	VALID_STRING_SOURCES = ['BLAST_KEGG_NAME','BioMart_HUGO']
	drugDict = {'sunitinib': {'source':['LUAD','SKCM']},
		'crizotinib':{'source':['LUAD','BRCA']},
		'sorafenib': {'source':['LUAD','SKCM','BRCA'],'target':['OV','LIHC','LUSC']},	
		'erlotinib': {'source':['LUAD','BRCA'],'target':['OV','LIHC','LUSC']}
			}

	annFile = "../data/raw/expression/CTRPv2/Cell_lines_annotations_20181226.txt"
	cclMetaFile = "../data/raw/expression/CTRPv2/v20.meta.per_cell_line.txt"
	
	compoundFile = "../data/raw/expression/CTRPv2/v20.meta.per_compound.txt"
	assayFile = "../data/raw/expression/CTRPv2/v20.meta.per_experiment.txt"
	responseFile = "../data/raw/expression/CTRPv2/v20.data.curves_post_qc.txt"

	tpmFile = "../data/raw/expression/CTRPv2/CCLE_RNAseq_rsem_genes_tpm_20180929.txt"
	countsFile = "../data/raw/expression/CTRPv2/CCLE_RNAseq_genes_counts_20180929.gct"
	featureFile = "../data/raw/expression/CTRPv2/CCLE_metabolomics_20190502.csv"

	cclData = fetch_ccl_annotations(annFile,cclMetaFile)
	
	
	compound, assay, responses = fetch_compound_info(
		compoundFile,
		assayFile,
		responseFile,
		list(drugDict.keys()))
	

	
	n_iters = 200
	
	print("\n\tReading TPM...")
	tpms = pd.read_csv(tpmFile,sep = "\t")

	print("\tReading Counts...")
	counts= pd.read_csv(countsFile,sep = "\t", skiprows = 2)

	print("\tReading Features..")
	feats= pd.read_csv(featureFile)


	kept_features = ['2-hydroxyglutarate','lactate', 'glutamine',
	 'asparagine', 'aspartate','tryptophan','methionine']
	
	feats = feats[['CCLE_ID']+kept_features]
	

	nameMap = counts[['Name','Description']]
	tpms.rename(columns = {'gene_id':'Name'},inplace = True)
	cclNames = list(tpms.columns)[2:]

	tpms = tpms.merge(nameMap,on="Name")
	tpms = tpms[["Description"]+cclNames]
	counts = counts.drop(columns=['Name'])
	counts.rename(columns = {'Description':'Gene'},inplace = True)
	tpms.rename(columns = {'Description':'Gene'},inplace = True)

	countsFull = counts.copy(deep = True)
	tpmsFull = tpms.copy(deep = True)


	measured_genes = list(counts['Gene'].values)

	STRING_aliases = pd.read_csv("../data/raw/networks/STRING/9606.protein.aliases.v12.0.txt",sep="\t")
	STRING_aliases = STRING_aliases[STRING_aliases['source'].isin(VALID_STRING_SOURCES)]
	STRING_links = pd.read_csv("../data/raw/networks/STRING/9606.protein.links.v12.0.txt",sep = ' ')
	
	good_pids, keep_genes, pid_gene_map = utils.find_unambiguous_nodes(STRING_aliases,measured_genes)
	

	print(len(keep_genes))
	sys.exit(1)

	
	count_genes_full = list(countsFull['Gene'].values)
	tpm_genes_full = list(tpmsFull['Gene'].values)

	
	countsFull = countsFull[countsFull['Gene'].isin(measured_genes)]

	countsFull.set_index('Gene',inplace = True)
	countsFull = countsFull.loc[~(countsFull==0).all(axis=1)]
	countsFull.reset_index(drop=False,inplace=True)
	countsFull = countsFull.groupby("Gene").mean().reset_index(drop=False)
	countsFull.drop_duplicates(subset=["Gene"],inplace=True)

	tpmsFull = tpms[tpms['Gene'].isin(measured_genes)]
	tpmsFull = tpmsFull[tpmsFull['Gene'].isin(list(countsFull['Gene'].values))]
	tpmsFull.set_index('Gene',inplace=True)
	tpmsFull = tpmsFull.loc[~(tpmsFull==0).all(axis=1)]
	tpmsFull.reset_index(drop=False,inplace=True)
	tpmsFull.drop_duplicates(subset=["Gene"],inplace=True)

	counts = counts[counts['Gene'].isin(keep_genes)]
	counts.set_index('Gene',inplace=True)
	counts = counts.loc[~(counts==0).all(axis=1)]
	counts.reset_index(drop=False,inplace=True)
	counts = counts.groupby("Gene").mean().reset_index(drop=False)
	counts.drop_duplicates(subset=["Gene"],inplace=True)
	
	
	tpms = tpms[tpms['Gene'].isin(keep_genes)]
	tpms = tpms[tpms['Gene'].isin(list(counts['Gene'].values))]
	tpms.set_index('Gene',inplace=True)
	tpms= tpms.loc[~(tpms==0).all(axis=1)]
	tpms.reset_index(drop=False,inplace=True)
	tpms.drop_duplicates(subset=["Gene"],inplace=True)

	counts = counts[counts['Gene'].isin(tpms['Gene'].values)]
	countsFull = countsFull[countsFull['Gene'].isin(tpmsFull['Gene'].values)]

	
	
	print("\tMaking Networks...")
	STRING_links = STRING_links[(STRING_links['protein1'].isin(good_pids)) & (STRING_links['protein2'].isin(good_pids))]
	
	os.makedirs("../data/processed/networks/ccle/",exist_ok = True)
	sparse, sparseW, dense, denseW, easy, easyW, tight, tightW=\
		utils.construct_networks(STRING_links,pid_gene_map, fname = "../data/processed/networks/ccle/stats_.txt")


	for network, name,stem in tqdm.tqdm(zip([sparse,sparseW,dense, denseW,easy, easyW,tight,tightW],
		['unweighted','weighted','unweighted','weighted','unweighted','weighted','unweighted','weighted'], 
		["sparse","sparse","dense","dense", "easy", "easy","tight", "tight"])):
		
		os.makedirs("../data/processed/networks/ccle/{s}".format(s=stem),exist_ok = True)
		with open("../data/processed/networks/ccle/{s}/{n}.pickle".format(s = stem,n =name),'wb') as ostream:
			pickle.dump(network,ostream)


	for drug in drugDict.keys():
		
		
		drugInfo = compound[compound['cpd_name']==drug]
		targets = drugInfo['gene_symbol_of_protein_target'].values[0].split(";")
		cpdID = drugInfo['master_cpd_id'].values[0]
		
		

		with open("../data/processed/genesets/{d}_targets.txt".format(d=drug),"w") as ostream:
			ostream.writelines([x+"\n" for x in targets])
		
		responseSubset = responses[responses['master_cpd_id']==cpdID].copy(deep=True)

		responseSubset.reset_index(inplace=True,drop=True)
		expData = responseSubset.merge(assay, on = "experiment_id")	
		
		
		cclResponse = expData[['master_ccl_id','area_under_curve']]

		cclResponse = cclResponse.groupby(['master_ccl_id']).mean().reset_index()
		
		
		vals = cclResponse['area_under_curve'].values
		q = np.quantile(vals, q = 1.0/3)


		cclResponse['Response'] = cclResponse['area_under_curve'].map(lambda x: 1 if x<=q else 0)
		cclResponse['Binned_Response'] = cclResponse['area_under_curve'].map(lambda x: "R" if x<=q else "NR")
		cclResponse['Thresh'] = q
		
		expData = cclResponse.merge(cclData, on = "master_ccl_id")

		
		sourceTissues = drugDict[drug]['source']
		
		sourceTissues =  sourceTissues + ['BROAD','PANCAN']
		# sourceTissues =  sourceTissues + ['PANCAN']

		for tissue in sourceTissues:
			
			os.makedirs("../data/processed/expression/ccle/{d}/{t}/".format(d=drug,t=tissue),exist_ok = True)
			
			if tissue != "PANCAN":
				tissueList = tissueDict[tissue]
			else:
				tissueList = drugDict[drug]['source']

			print("starting {d}\t{t}".format(d=drug,t=tissue))
			print(tissueList)

			tissueData = expData[expData['tcga_code'].isin(tissueList)]
			tissueCellLines = list(tissueData['CCLE_ID'].values)
			tissueCellLines = [x for x in tissueCellLines if x in counts.columns]
			
			colData = tissueData[['CCLE_ID','Binned_Response']]
			colData = colData[colData['CCLE_ID'].isin(tissueCellLines)]
			

			colData.rename(columns = {'CCLE_ID':'Run_ID'},inplace = True)
			
			colData.set_index('Run_ID',inplace = True)
			colData = colData.loc[tissueCellLines,:]
			colData.reset_index(inplace = True, drop = False, names = ['Run_ID'])

			countData = counts[['Gene']+tissueCellLines]

			countData.rename(columns = {'Gene':'Gene_ID'},inplace = True)
			

			colData.to_csv("./ccle_cols.csv",index = False)
			countData.to_csv("./ccle_counts.csv",index = False)


			countData.rename(columns = {'Gene_ID':'Gene'},inplace = True)
			

			tissueResponses = expData[['CCLE_ID','Response']]
			tissueResponses.rename(columns = {'CCLE_ID':'Run_ID'},inplace = True)
			tissueResponses = tissueResponses[tissueResponses['Run_ID'].isin(tissueCellLines)]
			tissueResponses.set_index('Run_ID',inplace = True)
			sampleOrder = list(tissueResponses.index)
			tissueResponses.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			tissueResponses.to_csv("../data/processed/expression/ccle/{d}/{t}/response.csv".format(d=drug,t=tissue),index=False)

			countsDataFull = countsFull[['Gene']+tissueCellLines]
			countsDataFull = countsDataFull.transpose()
			countsDataFull.reset_index(inplace = True)
			countsDataFull.columns = countsDataFull.iloc[0]

			countsDataFull.drop(index=0,axis=0,inplace=True)
			countsDataFull.rename(columns={'Gene':'Run_ID'},inplace = True)

			countsDataFull.set_index('Run_ID',inplace = True)
			countsDataFull = countsDataFull.loc[sampleOrder,:]
			countsDataFull.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			countsDataFull.to_csv("../data/processed/expression/ccle/{d}/{t}/counts_full.csv".format(d=drug,t=tissue),index = False)
			


			countData = countData.transpose()
			countData.reset_index(inplace = True)
			countData.columns = countData.iloc[0]
			countData.drop(index=0,axis=0,inplace=True)
			countData.rename(columns={'Gene':'Run_ID'},inplace = True)
			countData.set_index('Run_ID',inplace = True)
			countData = countData.loc[sampleOrder,:]
			countData.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			countData.to_csv("../data/processed/expression/ccle/{d}/{t}/counts.csv".format(d=drug,t=tissue),index = False)


			expressionFull = tpmsFull[['Gene']+tissueCellLines]
			expressionFull = expressionFull.transpose()
			expressionFull.reset_index(inplace = True)
			expressionFull.columns = expressionFull.iloc[0]
			expressionFull.drop(index=0,axis=0,inplace=True)
			expressionFull.rename(columns={'Gene':'Run_ID'},inplace = True)
			expressionFull.set_index('Run_ID',inplace = True)
			expressionFull = expressionFull.loc[sampleOrder,:]
			expressionFull.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			expressionFull.to_csv("../data/processed/expression/ccle/{d}/{t}/expression_full.csv".format(d=drug,t=tissue),index = False)


			expressionMeasurements = tpms[['Gene']+tissueCellLines]
			expressionMeasurements = expressionMeasurements.transpose()
			expressionMeasurements.reset_index(inplace = True)
			expressionMeasurements.columns = expressionMeasurements.iloc[0]
			expressionMeasurements.drop(index=0,axis=0,inplace=True)
			expressionMeasurements.rename(columns={'Gene':'Run_ID'},inplace = True)
			expressionMeasurements.set_index('Run_ID',inplace = True)
			expressionMeasurements = expressionMeasurements.loc[sampleOrder,:]
			expressionMeasurements.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			expressionMeasurements.to_csv("../data/processed/expression/ccle/{d}/{t}/expression.csv".format(d=drug,t=tissue),index = False)
			
			
			
			features = feats[feats['CCLE_ID'].isin(tissueCellLines)]
			features.rename(columns = {'CCLE_ID':'Run_ID'},inplace = True)
			issues = [x for x in sampleOrder if x not in list(features['Run_ID'])]
			print(issues)
			sampleOrder = [x for x in sampleOrder if x in list(features['Run_ID'])]
			features.set_index('Run_ID',inplace = True)
			features = features.loc[sampleOrder,:]
			features.reset_index(inplace = True, drop = False, names = ['Run_ID'])
			features.to_csv("../data/processed/expression/ccle/{d}/{t}/features.csv".format(d=drug,t=tissue),index = False)


			print("\n")
			print("Full Expression:\t{s}".format(s=expressionFull.shape))
			print("Full Count:\t{s}".format(s=countsDataFull.shape))
			print("Expression TPM:\t{s}".format(s=expressionMeasurements.shape))
			print("Expression Counts:\t{s}".format(s=countData.shape))
			print("Metab Features:\t{s}".format(s=features.shape))
			print("\n")
			
			

			if tissue != 'BROAD':
				cmd = "Rscript compute_DE_genes.R {drug} {tissue} {n_iters} ./ccle_counts.csv ./ccle_cols.csv".format(drug=drug, 
						tissue = tissue, 
						n_iters = n_iters)
			
				os.system(cmd)
		

			

	

	

	



if __name__ == '__main__':
	main()