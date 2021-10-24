# redologos
A **reductive Tricarboxilyc Acid Cycle** network evolutionary analysis.

[join our discord server](https://discord.gg/tHfetKqn)

>### Current dependencies
These are external dependencies for [redologos/bin](redologos/bin) scripts, links to their installing instructions, commands for installation,and used versions:
 | Software name (version used) 	| Terminal 	| Installation *debian based distros 	|
|-	|-	|-	|
| [tidyverse](https://www.tidyverse.org/) (1.3.1) 	| R (4.1.1) 	| `install.packages("tidyverse")` 	|
| [taxize](https://cran.r-project.org/web/packages/taxize/taxize.pdf) (0.9.99) 	| R (4.1.1) 	| `install.packages("taxize")` 	|
| [usethis](https://usethis.r-lib.org/) (2.0.1) 	| R (4.1.1) 	| `install.packages("usethis")` 	|
| [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) (4.2.8) 	| R (4.1.1) 	| `install.packages("seqinr")` 	|
| [myTAI](https://cran.r-project.org/web/packages/myTAI/index.html) (0.9.3) 	| R (4.1.1) 	| `install.packages("myTAI")` 	|
| [stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0) (1.4.0) 	| R (4.1.1) 	| `install.packages("stringr")` 	|
| [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html) (1.14.0) 	| R (4.1.1) 	| `install.packages("data.table")` 	|

## [bin/](bin/)
[Here](bin/) you'll find all scripts in the project. As a general rule, scripts in [bin/](bin/), read inputs from [data/](data/) and write outputs into [results/](results/). However external databases may be needed locally, and might be downloaded to independent locations.

Usage examples of scripts and more details with samples can be found in [bin/README.md](bin/README.md).

In brief, these are the current scripts in bin, it's ins-outs and tasks:

 | Script	| ins	| outs | tasks
|-	|-	|-	|-	|
| [Subset_by_pfam.R](bin/Subset_by_pfam.R) 	|[query.fasta](data/Proteins_fasta/All_queries.faa)<br />[subject.fasta](data/Proteins_fasta/All_Results_fastacmd3.faa)<br />path to Pfamm-A.hmm | [Subset.fasta](results/subset_by_pfam/All_Results_fastacmd3_subset_by_All_queries_pfams.fasta)<br />[Subset.csv](results/subset_by_pfam/All_Results_fastacmd3_subset_by_All_queries_pfams.csv) | Inkvokes **hmmscan** to list all Pfam IDs in query.fasta, then retrieves all sequences in subject.fasta containing at least one of listed Pfam IDs
