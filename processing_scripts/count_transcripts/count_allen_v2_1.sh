# generate count matrices 221116
# lets quantify some v2 allen data

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv2 \
-o ./allen_B05/ \
-t 20 -m 50G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/allen_B05/L8TX_171026_01_B05_R1.fastq.gz   \
../datasets/allen_B05/L8TX_171026_01_B05_R2.fastq.gz 
