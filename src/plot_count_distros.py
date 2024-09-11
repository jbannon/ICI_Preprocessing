from scipy.stats import kurtosis,mode, skew, beta
import numpy as np
import os 
import sys
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import diptest
from collections import defaultdict



goodThresh = 2.0/3.0
ProbThresh = 0.90
DRUG_TISSUE_MAP = {
	"Atezo":["KIRC","BLCA"],
	"Pembro":["SKCM","STAD"],
    "Nivo":["SKCM","KIRC"], 
    "Ipi":["SKCM"], 
    "Ipi+Pembro":["SKCM"],
    'sorafenib':["BRCA","LUAD","SKCM"],
    'erlotinib':['LUAD','BRCA'],
    'sunitinib': ['LUAD','SKCM'],
	'crizotinib':['LUAD','BRCA']}
DRUG_TISSUE_MAP = {

    "Nivo":["SKCM","KIRC"], 
 }

for drug in DRUG_TISSUE_MAP.keys():
	
	for tissue in DRUG_TISSUE_MAP[drug]:
	
		de_file = "../data/processed/genesets/{d}_{t}_DE.csv".format(d=drug,t=tissue)	
		

		df = pd.read_csv(de_file)
	
		
		#["tight","narrow","very.tight"]:
		

		for fdr in list(pd.unique(df["Thresh.Name"])):
			
			print("************************")
			print("\t{d} - {t} - {f}".format(d=drug,t=tissue,f=fdr))
			print("************************")

			
			res = "./count_figs/full/{d}/{t}/{f}/".format(d=drug,t=tissue,f=fdr)
			os.makedirs(res,exist_ok=True)
			# print(df)
			temp = df[df['Thresh.Name']==fdr]
			
			val = temp['Thresh.Value'].values[0]
			
			ax = sns.histplot(data=temp, x= 'Count')#,kde=True)
			
			
			
			# temp = temp[temp['Count']>=0]
			N = len(pd.unique(temp['Gene']))
			temp["Hits"] = temp['Count']
			temp['Misses'] = 200-temp['Count']
			# temp = temp[["Gene","Count"]]
			temp['Proportion'] = temp["Count"]/200
			# print(temp["Proportion"])
			
			# print(np.amax(temp["Proportion"].values))
			a, b,loc, scale = beta.fit(temp[temp["Proportion"]<1.0]["Proportion"].values,fscale = 1, floc = 0)
			# print("tong")
			temp["EB_avg"] = (temp["Hits"]+a)/(200+a+b)
			print(temp)
			continue
			EB_point = np.amin(temp[temp["EB_avg"]>=goodThresh]["Count"].values)
			
			
			prob_estimates = defaultdict(list)
			for idx, row in temp.iterrows():
				
				a1 = a + row["Hits"]
				b1 = b + row["Misses"]
				thisBeta = beta(a=a1,b=b1)
				thisProb = thisBeta.sf(goodThresh)
				prob_estimates['Gene'].append(row['Gene'])
				prob_estimates['EB_Prob'].append(thisProb)
				
			
			pb = pd.DataFrame(prob_estimates)

			EB = pb[pb['EB_Prob']>=ProbThresh]
			
			EB = temp[temp["Gene"].isin(EB["Gene"].values)]
			try:
				EB_cut = np.amin(EB["Count"].values)
			except: 
				EB_cut = 200
			min_cut = min(EB_point,EB_cut)

			
			num_genes = len(pd.unique(temp['Gene']))
			ax = sns.histplot(data=temp, x= 'Count')#,kde=True)
			ax.axvline(EB_point, color = "red",label = "EB-Point")
			ax.axvline(EB_cut, color = "magenta",label = "EB-CDF")
			
			ax.set(title = "{d} {t}. Pvalue thresh: {pval}, Total Genes: {n}".format(d=drug,t=tissue,n=N,pval=val))
			plt.legend()
			
			plt.tight_layout()
			plt.savefig("{r}{d}_{t}_{f}.png".format(r=res,d=drug,t=tissue,f=fdr))
			plt.close()
			
			temp = temp[temp['Count'] >= min_cut]
			num_genes = len(pd.unique(temp['Gene']))
			ax = sns.histplot(data = temp, x='Count')#,kde=True)
			ax.axvline(EB_point, color = "red",label = "EB-Point")
			ax.axvline(EB_cut, color = "magenta",label = "EB-CDF")
			ax.set(title = "{d} {t} {cut}\n {n} Total Genes, Trimmed".format(d=drug,t=tissue, cut=val, n = num_genes))
			plt.legend()
			plt.tight_layout()
			plt.savefig("{r}{d}_{t}_{f}_trimmed.png".format(r=res,d=drug,t=tissue,f=fdr))
			plt.close()
			
			
	




