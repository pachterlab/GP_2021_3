#!/bin/bash

# download  10x nuclei-matched full cell dataset: 210622
mkdir -p brain_10x_5k_fastqs
wget -P ./brain_10x_5k_fastqs https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_Neurons_5K_Multiplex/SC3_v3_NextGem_DI_Neurons_5K_Multiplex_fastqs.tar
tar -xvf ./brain_10x_5k_fastqs/*.tar
