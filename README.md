# SIDR2.0 (under construction)

"Sequence Identification with Decision tRees"

SIDRv2 is an expanded version of the original SIDR method: a sequence identification/decontamination tool that classifies unknown genetic material collected from environmental samples.

The majority of SIDR2.0 is written in C to decrease execution time of data parsing, while machine learning models remain programmed in Python using Scikit-learn. SIDR2.0 utilizes the Python3 developer package to interface between C and Python, therefore allowing increased runtime efficiency. We're also adding an extra model for expectation maximization of kmer distribution data to increase robustness of decision tree classification and eliminate classification errors.

Original SIDR method: https://doi.org/10.1186/s12859-017-1941-0

Kmer distribution algorithm: https://doi.org/10.1093/bib/bbaa063

Sequencing data backend API: https://github.com/samtools/htslib
