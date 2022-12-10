## gRNA Designer

### Description

This script designated for searching and evaluating the candidate gRNAs. Based on a simple motivation, only two files of target sequence and genome sequence are needed to automatically design gRNAs.

### Design pipeline
We use [`cas-offinder`](https://github.com/snugel/cas-offinder) as search engine to search each possible candidate in **Target sequence**, and than we evaluate each of them in **Score sequence**.


### Usage
Modify and Run.

```python
design_path = {Path of Target sequence, .fasta file}
score_path = {Path of Score sequence, .fasta file}
Cas_species = 'SpCas9'  # Which species of Cas protein used in project.
task_name = None
work_dir = None

DNAbuldge = 1
RNAbuldge = 2
mismatch = 5
risk_match = 3

------------
```
