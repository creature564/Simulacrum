# SIML (under construction)

"Sequence Identification using Machine Learning"

SIML is an expanded version of the original SIDR method: a sequence identification/decontamination tool that classifies unknown genetic material collected from environmental samples.

The majority of SIML is written in C to decrease execution time of data parsing, while machine learning models remain programmed in Python using Scikit-learn. SIML utilizes the Python3 developer package to interface between C and Python, therefore allowing increased runtime efficiency. We're also adding an extra model for expectation maximization of kmer distribution data to increase robustness of decision tree classification and decrease classification errors.

Original SIDR method: https://doi.org/10.1186/s12859-017-1941-0

Kmer distribution algorithm: https://doi.org/10.1093/bib/bbaa063

Sequencing data backend API: https://github.com/samtools/htslib
