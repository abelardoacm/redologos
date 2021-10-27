# Bin navigation details

In this file you'll find details of the individual steps, what the scripts are doing and their inputs and outputs so you know how each folder is generated.

## **Subset_by_pfam.R**
Task: Inkove hmmscan twice (if needed) to scan two fasta files searching for Pfam domains. When Pfam IDs for query and subject are obtained, then a subset of subject is built, keeping only those sequences with Pfam IDs listed for query sequences.

>#### Sample input:
- **query** : [All_queries.faa](/data/Proteins_fasta/All_queries.faa)
- **subject** : [All_Results_fastacmd3.faa](/data/Proteins_fasta/All_Results_fastacmd3.faa)
- **path to Pfamm-A.hmm** : */home/abelardo/get_homologues/db/Pfam-A.hmm*

>#### Command
```
$ Rscript Subset_by_pfam.R -q ../data/Proteins_fasta/All_queries.faa -s ../data/Proteins_fasta/All_Results_fastacmd3.faa -P /home/abelardo/get_homologues/db/Pfam-A.hmm
```
>#### stdout
```
Previous hmmscan found for query!
Previous hmmscan found for subject!

Number of proteins in query "All_queries" = 7	|	Number of Pfams = 13
Number of proteins in subject "All_Results_fastacmd3" = 439	|	Number of Pfams = 36

Subject counts of sequences with Pfams matching query:
     23 "ATP_citrate_lyase_citrate-binding"
     29 "ATP-grasp_domain"
     26 "Citrate_synthase
     29 "CoA_binding_domain"
     27 "CoA-ligase"
     62 "Pyruvate_ferredoxin/flavodoxin_oxidoreductase"
     50 "Pyruvate_ferredoxin_oxidoreductase_beta_subunit_C_terminal"
     68 "Pyruvate:ferredoxin_oxidoreductase_core_domain_II"
     98 "Pyruvate_flavodoxin/ferredoxin_oxidoreductase
     14 "Succinyl-CoA_ligase_like_flavodoxin_domain"
     70 "Thiamine_pyrophosphate_enzyme
```
>#### Sample output:
- **subset of subject sequences in fasta**: [All_Results_fastacmd3_subset_by_All_queries_pfams.fasta](results/subset_by_pfam/All_Results_fastacmd3_subset_by_All_queries_pfams.fasta)
- **subject hmmscan out**: [All_Results_fastacmd3_subset_by_All_queries_pfams.csv](results/subset_by_pfam/All_Results_fastacmd3_subset_by_All_queries_pfams.csv)

>#### Arguments:
- Needed
    - **-q** \<path to query.fasta\>
    - **-s** \<path to subject.fasta\>
    - **-P** \<path to Pfam\>
- Optional
    - **--stdout** prints Pfam fasta subset in stdout
    - **-o** optional name for output instead of default (subject_subset_by_query_pfams.csv)
