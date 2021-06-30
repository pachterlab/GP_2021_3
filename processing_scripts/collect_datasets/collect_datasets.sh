#!/bin/bash

# download 10x datasets: 210517

curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar
tar -xvf pbmc_10k_v3_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar
tar -xvf neuron_1k_v3_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_fastqs.tar
tar -xvf neuron_10k_v3_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_1k_v3/heart_1k_v3_fastqs.tar
tar -xvf heart_1k_v3_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/heart_10k_v3/heart_10k_v3_fastqs.tar
tar -xvf heart_10k_v3_fastqs.tar
