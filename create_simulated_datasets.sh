#!/bin/bash
# S.agalactiae
declare -a cluster100=("ERR3468293" "SRR7283372" "SRR2062125" "SRR7281691" "ERR1910678" )
declare -a cluster30=("SRR5063572" "ERR1741401" "SRR5061557" "SRR3320830" "ERR1659686" )
declare -a cluster20=("ERR1741950" "SRR7282364" "SRR6050836" "ERR1672714" "ERR3467983" )

for val in ${cluster100[@]};do
    seqtk sample -s100 ${val}_1.fastq 100000 > ${val}_1.100.fastq
    seqtk sample -s100 ${val}_2.fastq 100000 > ${val}_2.100.fastq
done

for val in ${cluster30[@]};do
    seqtk sample -s100 ${val}_1.fastq 30000 > ${val}_1.30.fastq
    seqtk sample -s100 ${val}_2.fastq 30000 > ${val}_2.30.fastq
done

for val in ${cluster20[@]};do
    seqtk sample -s100 ${val}_1.fastq 20000 > ${val}_1.20.fastq
    seqtk sample -s100 ${val}_2.fastq 20000 > ${val}_2.20.fastq
done
