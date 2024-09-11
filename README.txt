This is a repository to process the CRI and GDSC data subsets. 

Usage:

to prepare datasets proper:

python3 prepare_cri.py -config configs/log_tpm.yaml 
python3 prepare_cri.py -config configs/tpm.yaml
python3 prepare_GDSC.py -config configs/GDSC.yaml

for genesets:

python3 prepare_genesets.py -config configs/genesets.yaml