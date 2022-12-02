#!/bin/bash

# generate count matrices 221115

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv2 \
-o ./neuron_1k_v2/ \
-t 16 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/neuron_1k_v2_fastqs/neuron_1k_v2_S1_L001_R1_001.fastq.gz \
../datasets/neuron_1k_v2_fastqs/neuron_1k_v2_S1_L001_R2_001.fastq.gz \
../datasets/neuron_1k_v2_fastqs/neuron_1k_v2_S1_L002_R1_001.fastq.gz \
../datasets/neuron_1k_v2_fastqs/neuron_1k_v2_S1_L002_R2_001.fastq.gz
