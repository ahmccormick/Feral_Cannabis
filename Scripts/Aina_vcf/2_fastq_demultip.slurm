#!/bin/bash
#SBATCH --job-name=demultip
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH -e /home/ahmccorm/%x-%A_%a.err
#SBATCH -o /home/ahmccorm/%x-%A_%a.out
#SBATCH --mail-user=ahmccorm@hawaii.edu
#SBATCH --mail-type=ALL

module load bio/Stacks/2.62-foss-2022a

grep -v '^#' /home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3/GBS_barcodes_clean.txt \
  > /home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3/GBS_barcodes_clean.nohdr.txt


IN=/home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3/adaptor_trimmed_fastq
OUT=/home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3/demultiplexed
BAR=/home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3/GBS_barcodes_clean.nohdr.txt

mkdir -p "$OUT"

process_radtags \
  -P -p "$IN" \
  -o "$OUT" \
  -b "$BAR" \
  -e apeKI \
  --inline_null \
  -i gzfastq \
  -c -q -r
