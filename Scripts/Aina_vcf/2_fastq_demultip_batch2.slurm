#!/bin/bash
#SBATCH --job-name=demultip
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH -e /home/ahmccorm/%x-%A_%a.err
#SBATCH -o /home/ahmccorm/%x-%A_%a.out
#SBATCH --mail-user=ahmccorm@hawaii.edu
#SBATCH --mail-type=ALL

module load bio/Stacks/2.62-foss-2022a

BASE=/home/ahmccorm/koa_scratch/230411_AHY33JDSX5
BAR=${BASE}/GBS_barcodes_clean.nohdr.txt   # two columns: <barcode>\t<well> (A01..H12)

for P in 1 2 3 4; do
  awk -v P=$P 'BEGIN{FS=OFS="\t"} {print $1, "Plate-" P "-" $2}' \
    "$BAR" > ${BASE}/barcodes_plate${P}.txt
done

IN=${BASE}/adaptor_trimmed_fastq
OUT=${BASE}/demultiplexed
mkdir -p "$OUT"

for P in 1 2 3 4; do
  R1=$(ls ${IN}/Plate-${P}_*_R1_trimmed.fastq.gz)
  R2=${R1/_R1_/_R2_}
  process_radtags \
    -P \
    -1 "$R1" \
    -2 "$R2" \
    -o "$OUT" \
    -b ${BASE}/barcodes_plate${P}.txt \
    -e apeKI \
    --inline_null \
    -i gzfastq \
    -c -q -r
done
