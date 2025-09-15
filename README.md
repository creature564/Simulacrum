# SIML -- Sequence Identification with Machine Learning

"Sequence Identification using Machine Learning"

SIML is an expanded version of the original SIDR method: a sequence identification/decontamination tool that classifies unknown genetic material collected from environmental samples.

This program is a demo version of the C backend code I used to achieve massive processing speed boosts on genomic data, as compared to standard Python libraries. Data is processed and cleaned in C, then ported to a Python environment for machine learning. The machine learning routines are not shown here, as they were directly rewritten from the original SIDR method. The major modification was writing the manual ETL for data processing and transmission across environments.

Original SIDR method: https://doi.org/10.1186/s12859-017-1941-0

Kmer distribution algorithm: https://doi.org/10.1093/bib/bbaa063

Sequencing data backend API: https://github.com/samtools/htslib
