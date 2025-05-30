'''
This Script was designed for search candidates of gRNA in Genome.

Cas_species list: SpCas9, AsCas12f, Un1Cas12f, FnCas12a

@author: CHU M. Pan
@Email: pan_chu@outlook.com
'''
# %%
import os
import json
import subprocess
import pandas as pd
import numpy as np
import argparse
import sys
CAS_OFFINDER_PATH = r'./cas-offinder'
TARGET_PATTERN_PATH= r'./targetRNAPattern.json'
# %% load default parameters
design_path = r'./Inputs/spCas9_pit.fa'
score_path = r'./Inputs/NCM3722.fasta'
Cas_species = 'SpCas9' #'SpCas9' 'AsCas12f' 'Un1Cas12f'
task_name = None  # if None, default is the input file name.
work_dir = None  # if not specify work directory, default set is the current working path and creating a new fold with name of the task name

DNAbuldge = 1
RNAbuldge = 2
mismatch = 5
risk_match = 3
# read input parameters
# -i <design_path> the file path of the designed gRNA, it's fasta format
# -s <score_path> the file path of the scoring sequences, it's fasta format
# -c <Cas_species> the Cas species to be used, default is SpCas9, we can use AsCas12f, Un1Cas12f, FnCas12, that is defined in the targetRNAPattern.json
# -t <task_name> the name of the task
# -w <work_dir> the working directory for output files
# -h or --help for printing help information
# help information: This script was designed for search candidates of gRNA in Genome.
# usage: python gRNADesigner.py -i <design_path> -s <score_path> -c <Cas_species> -t <task_name> -w <work_dir>
parser = argparse.ArgumentParser(description='This script was designed for search candidates of gRNA in Genome.')
parser.add_argument('-i', '--input_file', type=str, default=design_path,
                    help='the file path of the designed gRNA, it\'s fasta format')
parser.add_argument('-s', '--score_path', type=str, default=score_path,
                    help='the file path of the scoring sequences, it\'s fasta format')
parser.add_argument('-c', '--Cas_species', type=str, default=Cas_species,
                    help='the Cas species to be used, default is SpCas9, we can use AsCas12f, Un1Cas12f, FnCas12, that is defined in the targetRNAPattern.json')
parser.add_argument('-t', '--task_name', type=str, default=task_name,
                    help='the name of the task, if None, default is the input file name.')
parser.add_argument('-w', '--work_dir', type=str, default=work_dir,
                    help='the working directory for output files, if not specify work directory, default set is the current working path and creating a new fold with name of the task name')
args = parser.parse_args()
# 1. read the input parameters
if args.input_file is not None:
    design_path = args.input_file
if args.score_path is not None:
    score_path = args.score_path
if args.Cas_species is not None:
    Cas_species = args.Cas_species
if args.task_name is not None:
    task_name = args.task_name
if args.work_dir is not None:
    work_dir = args.work_dir
# 2. check the input file
if not os.path.exists(design_path):
    print(f'Error: The design file {design_path} does not exist.')
    sys.exit(1)
# 3. user check the input file path and parameters for the script
# 3.a print info
print(f'Input design file: {design_path}')
print(f'Input score file: {score_path}')
print(f'Cas species: {Cas_species}')
print(f'Task name: {task_name}')
while True:

    
    user_input = input("Please confirm the design path and parameters (y/n): ")
    if user_input.lower() == 'y':
        break
    elif user_input.lower() == 'n':
        print("Exiting the script. Please check your inputs.")
        sys.exit(1)

# %%

if task_name is None:
    design_name = os.path.basename(design_path)
    design_name = design_name.split('.')[:-1]
    if len(design_name) > 1:
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
casOfFinderPath = CAS_OFFINDER_PATH
pam_dict = json.load(open(TARGET_PATTERN_PATH))
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
search_cmd = f'{casOfFinderPath} {search_file} C0 {search_outPutFile}'
print(f'[{task_name}] -> {search_cmd}')
subprocess.run(search_cmd, shell=True)
# %% 
# 1. Read all searched candidates
targets_candidates = pd.read_csv(search_outPutFile, sep='\t', header=None)
targets_seq = targets_candidates[3]

targets_filter = [True if isinstance(seq, str) else False for seq in targets_seq]
targets_seq = targets_candidates[3][targets_filter]
targets_strand = targets_candidates[4][targets_filter]
targets_loc = targets_candidates[2][targets_filter]
targets_seq_name = targets_candidates[1][targets_filter]  # the sequence targeted sequence name
pam_length = len(pam_dict[Cas_species]['PAM'])
pattern_length = len(pam_dict[Cas_species]['Pattern'])

pam_searchPattern = 'N' * pam_length
if pam_dict[Cas_species]['PAM_loc'] == 5:
    targets_seq = [seq[pam_length:pam_length + pattern_length] for seq in targets_seq]
    score_seqs = [pam_searchPattern + seq for seq in targets_seq]
else:
    targets_seq = [seq[:pattern_length] for seq in targets_seq]
    score_seqs = [seq + pam_searchPattern for seq in targets_seq]

# evaluate all search candidates
score_file = score_path + '\n' + \
             f'{findPattern} {DNAbuldge} {RNAbuldge}\n' + \
             '\n'.join([f'{seq} {mismatch}' for seq in score_seqs])

score_file_path = os.path.join(work_dir, f'{task_name}_score_gRNA.input')
with open(score_file_path, 'a') as file:
    file.write(score_file)

score_outPutFile = os.path.join(work_dir, f'{task_name}_score_gRNA.output')
score_cmd = f'{casOfFinderPath} {score_file_path} C0 {score_outPutFile}'
print(f'[{task_name}] -> {score_cmd}')
result = subprocess.run(score_cmd, shell=True)
# wait the score_cmd finish
# Check the result
if result.returncode == 0:
    print("Command executed successfully")
else:
    print("Command failed")

# %% Statistic for each designed sequence
score_ret = pd.read_csv(score_outPutFile, sep='\t', header=None)
GCs = []
matcheds_num = []
identicals_num = []
genome_loc = []
risks_num = []
for seqI in range(len(targets_seq)):
    target = targets_seq[seqI]
    # target_name = targets_seq_name[seqI]
    GC_ratio = np.sum([True if alpha in 'GC' else False for alpha in target]) / len(target) * 100.
    GCs.append(GC_ratio)
    score_seq = score_seqs[seqI]
    matched_seqs = score_ret[score_ret[0] == score_seq]
    matched_num = len(matched_seqs)
    matcheds_num.append(matched_num)
    identical_num = np.sum(matched_seqs[5] == 0)
    identicals_num.append(identical_num)
    if identical_num > 0:
        identical_loc = matched_seqs[matched_seqs[5] == 0][2].tolist()
        identical_loc = [str(item) for item in identical_loc]
        if identical_num == 1:
            genome_loc.append(identical_loc[0])
        else:
            genome_loc.append('; '.join(identical_loc))
    else:
        genome_loc.append(' ')
    risk_num = np.sum(np.logical_and(matched_seqs[5] <= risk_match, matched_seqs[5] > 0))
    risks_num.append(risk_num)

# # compare the length of all sequence
# print(f'Target: {len(targets_seq)}')
# print(f'TargetName: {len(targets_seq_name)}')
# print(f'Strand: {len(targets_strand)}')
# print(f'GC: {len(GCs)}')
# print(f'LocationinDesign: {len(targets_loc)}')
# print(f'Identical: {len(identicals_num)}')
# print(f'LocationinScore: {len(genome_loc)}')
# print(f'Match: {len(matcheds_num)}')
# print(f'Risk: {len(risks_num)}')


data = pd.DataFrame(data=dict(Target=targets_seq,
                              TargetName = targets_seq_name,
                              Strand=targets_strand,
                              GC=GCs,
                              LocationinDesign=targets_loc,
                              Identical=identicals_num,
                              LocationinScore=genome_loc,
                              Match=matcheds_num,
                              Risk=risks_num))

data.to_excel(os.path.join(work_dir, f'{task_name}_statistics.xlsx'))
# %%
12369874


