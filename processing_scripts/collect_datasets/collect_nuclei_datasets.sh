#!/bin/bash

# download  10x nuclei dataset: 210616
mkdir -p brain_nuc_10x_5k_fastqs
wget -P ./brain_nuc_10x_5k_fastqs https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex/SC3_v3_NextGem_DI_Nuclei_5K_Multiplex_fastqs.tar
tar -xvf ./brain_nuc_10x_5k_fastqs/*.tar
