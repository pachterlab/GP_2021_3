
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7388133.svg)](https://doi.org/10.5281/zenodo.7388133)

# Analysis scripts for "Length Biases in Single-Cell RNA Sequencing of pre-mRNA"

## Overview

This repository contains all of the quantification, inference, analysis, and visualization code for "Length Biases in Single-Cell RNA Sequencing of pre-mRNA". We investigate a variety of 10x v2 and v3 datasets, and find a pervasive length bias in the unspliced counts. Upon fitting a two-stage model of bursty transcription, splicing, and degradation, we discover non-physical trends in parameter values. We provide an explanatory model by proposing that the sequencing of unspliced RNA exhibits a length-dependent bias. Finally, we use the results from this model to investigate the relationships between paired datasets, including the characterization of differences between v2 and v3 technologies and the investigation of marker-like genes that show considerable distributional differences between cell types.

All analyses use [Monod](https://github.com/pachterlab/monod) version 0.2.5.0. For [documentation](https://monod-examples.readthedocs.io/en/latest/) and [examples](https://github.com/pachterlab/monod_examples), see the corresponding pages.

All raw data are available on [Zenodo](https://zenodo.org/record/7388133).

## Repository contents 

* `polyA_ref/`: gene lengths and numbers of poly(A) stretches for the human and mouse genomes.
* `processing_scripts/`: scripts used to perform preliminary data analysis.
  * `collect_datasets/`: scripts used to download some of the analyzed datasets. Note that most datasets were downloaded from the command line. The data locations are reported in the article. 
  * `count_transcripts/`: scripts used to pseudoalign all datasets.
  * `make_references/`: scripts used to generate genome references, and notebooks used to extract gene lengths and poly(A) stretch counts. The raw `FASTA` files used for the latter analysis can be obtained from Ensembl BioMart by querying gene name, gene start, and gene end properties for Attributes -> Sequences -> Unspliced (Gene).
* `notebooks/`: notebooks used to perform data analysis.
  *  `noisefree_10x_pbmc_fit`, `noisefree_10x_heart_fit`, `noisefree_10x_neuron_fit`, and `noisefree_allen_fit`: fits to entire datasets with a model that does not include technical noise.
  *  `10x_pbmc_fit`, `10x_heart_fit`, `10x_neuron_fit`, and `allen_fit`: grid fits to entire datasets using a length-biased model of technical noise.
  *  `celltypes_10x_pbmc_fit` and `celltypes_allen_fit`: fits to cell type populations using a length-biased model of technical noise.
  *  `MoM_10x_pbmc_fit`: a benchmark for the method of moments initialization used for inference.
  *  `nolengths_10x_pbmc_fit`: a fit using a length-independent model of technical noise.
  *  `figure_generation`: analysis and visualization of results produced by all of the other notebooks.
  *  `figs/`: figures output by `figure_generation.ipynb`.

## Software and hardware requirements

Data analysis was performed using up to 60 CPUs (3.7GHz each) on a dedicated server. The analyses used the following software:

```
python 3.9.1
scipy 1.9.1
numpy 1.23.2
matplotlib 3.6.0
loompy 3.0.6
numdifftools 0.9.40
anndata 0.8.0
kb-python 0.26.0
monod 0.2.5.0
```
