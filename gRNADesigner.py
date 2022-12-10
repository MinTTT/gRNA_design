'''
This Script was designed for search candidates of gRNA in Genome.

Cas_species list: SpCas9, AsCas12f, Un1Cas12f, FnCas12a

@author: CHU M. Pan
@Email: pan_chu@outlook.com
'''
#%%
import os
import json
import subprocess
import pandas as pd
import numpy as np
#%%
design_path = r'./Inputs/CPP442_promoter.fa'
score_path = r'./NCM3722.fasta'
Cas_species = 'SpCas9'
task_name = None
work_dir = None

DNAbuldge = 1
RNAbuldge = 2
mismatch = 5
risk_match = 3
#%%

if task_name is None:
    design_name = os.path.basename(design_path)
    design_name = design_name.split('.')[:-1]
    if len(design_name) >1:
        design_name = '.'.join(design_name)
    else:
        design_name = design_name[0]
    task_name = f'{Cas_species}_{design_name}'

if work_dir is None:
    try:
        work_dir = f'./{task_name}'
        os.makedirs(work_dir)
    except FileExistsError:
        print('Cannot create a file when that file already exists')
        exit(1)

# %%
pam_dict = json.load(open(r'./targetRNAPattern.json'))
PAM_seq = pam_dict[Cas_species]['PAM']
if pam_dict[Cas_species]['PAM_loc'] == 5:
    findPattern = PAM_seq + pam_dict[Cas_species]['Pattern']
else:
    findPattern = pam_dict[Cas_species]['Pattern'] + PAM_seq


find_file = design_path + '\n' + \
        f'{findPattern} 0 0\n' + \
        f'{findPattern} 0'

search_file = os.path.join(work_dir, f'{task_name}_search_gRNA.input')
with open(search_file, 'a') as file:
    file.write(find_file)

search_outPutFile = os.path.join(work_dir, f'{task_name}_search_gRNA.output')
search_cmd = f'./cas-offinder {search_file} G0 {search_outPutFile}'
print(f'[{task_name}] -> {search_cmd}')
subprocess.run(search_cmd, shell=True)
# %%
targets_candidates = pd.read_csv(search_outPutFile, sep='\t', header=None)
targets_seq = targets_candidates[3]
targets_filter = [True if isinstance(seq, str) else False for seq in targets_seq]
targets_seq = targets_candidates[3][targets_filter]
targets_strand = targets_candidates[4][targets_filter]
targets_loc = targets_candidates[2][targets_filter]
pam_length = len(pam_dict[Cas_species]['PAM'])
pattern_length = len(pam_dict[Cas_species]['Pattern'])

pam_searchPattern = 'N'*pam_length
if pam_dict[Cas_species]['PAM_loc'] == 5:
    targets_seq = [seq[pam_length:pam_length+pattern_length] for seq in targets_seq]
    score_seqs = [pam_searchPattern + seq for seq in targets_seq]
else:
    targets_seq = [seq[:pattern_length] for seq in targets_seq]
    score_seqs = [seq + pam_searchPattern for seq in targets_seq]

score_file = score_path + '\n' + \
        f'{findPattern} {DNAbuldge} {RNAbuldge}\n' + \
        '\n'.join([f'{seq} {mismatch}' for seq in score_seqs])

score_file_path = os.path.join(work_dir, f'{task_name}_score_gRNA.input')
with open(score_file_path, 'a') as file:
    file.write(score_file)

score_outPutFile = os.path.join(work_dir, f'{task_name}_score_gRNA.output')
score_cmd = f'./cas-offinder {score_file_path} G0 {score_outPutFile}'
print(f'[{task_name}] -> {score_cmd}')
subprocess.run(score_cmd, shell=True)
# %% Statistic for each designed sequence

score_ret = pd.read_csv(score_outPutFile, sep='\t', header=None)
GCs = []
matcheds_num = []
identicals_num = []
genome_loc = []
risks_num = []
for seqI in range(len(targets_seq)):
    target = targets_seq[seqI]
    GC = np.sum([True if alpha in 'GC' else False for alpha in target]) / len(target) * 100.
    GCs.append(GC)
    score_seq = score_seqs[seqI]
    matched_seqs = score_ret[score_ret[0] == score_seq]
    matched_num = len(matched_seqs)
    matcheds_num.append(matched_num)
    identical_num = np.sum(matched_seqs[5] ==0)
    identicals_num.append(identical_num)
    if identical_num > 0:
        identical_loc = matched_seqs[matched_seqs[5]==0][2].tolist()
        identical_loc = [str(item) for item in identical_loc]
        if identical_num == 1:
            genome_loc.append(identical_loc[0])
        else:
            genome_loc.append('; '.join(identical_loc))
    else:
        genome_loc.append(' ')
    risk_num = np.sum(np.logical_and(matched_seqs[5] <= risk_match, matched_seqs[5] > 0))
    risks_num.append(risk_num)

data = pd.DataFrame(data=dict(Target=targets_seq, 
                                Strand=targets_strand,
                                GC=GCs,
                                LocationinDesign=targets_loc,
                                Identical=identicals_num,
                                LocationinScore=genome_loc,
                                Match=matcheds_num,
                                Risk=risks_num))

data.to_excel(os.path.join(work_dir, f'{task_name}_statistics.xlsx'))
# %%
