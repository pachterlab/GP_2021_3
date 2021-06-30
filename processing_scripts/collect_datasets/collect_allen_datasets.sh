#!/bin/bash

# download allen datasets: 210604
#mkdir -p allen_B01 allen_C01 allen_A08 allen_B08
wget -P ./allen_B01 https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/raw/MOp/L8TX_181211_01_B01_S01_L003.fastq.tar
#tar -xvf ~/L8TX_181211_01_B01_S01_L003.fastq.tar

wget  https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/raw/MOp/L8TX_181211_01_C01_S01_L003.fastq.tar
#tar -xvf ~/L8TX_181211_01_C01_S01_L003.fastq.tar

wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/raw/MOp/L8TX_190430_01_A08_S01_L003.fastq.tar
#tar -xvf ~/L8TX_190430_01_A08_S01_L003.fastq.tar

wget https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/raw/MOp/L8TX_190430_01_B08_S01_L003.fastq.tar
#tar -xvf ~/L8TX_190430_01_B08_S01_L003.fastq.tar
