#!/bin/bash
#SBATCH --job-name=adaptor_trim
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
# NOTE: do NOT put --array here; pass it on the sbatch command line


module load bio/cutadapt/4.2-GCCcore-11.3.0

INDIR="/home/ahmccorm/koa_scratch/20240529_LH00283_0104_B22K7G3LT3"
OUTDIR="${INDIR}/adaptor_trimmed_fastq"
mkdir -p "$OUTDIR"

# build list of R1 files
#mapfile -t R1_FILES < <(ls "${INDIR}"/*_R1_001.fastq.gz)
mapfile -t R1_FILES < <(ls -1 "${INDIR}"/*_R1_001.fastq.gz | sort)


R1_FILE="${R1_FILES[$SLURM_ARRAY_TASK_ID]}"
R2_FILE="${R1_FILE/_R1_/_R2_}"

BASENAME="$(basename "$R1_FILE" _R1_001.fastq.gz)"
OUT_R1="${OUTDIR}/${BASENAME}_R1_trimmed.fastq.gz"
OUT_R2="${OUTDIR}/${BASENAME}_R2_trimmed.fastq.gz"
REPORT="${OUTDIR}/${BASENAME}_cutadapt.txt"

# skip if already done
if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
  echo "[$BASENAME] outputs exist, skipping."
  exit 0
fi

cutadapt \
  -j ${SLURM_CPUS_PER_TASK} \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -m 30 \
  -o "$OUT_R1" -p "$OUT_R2" \
  "$R1_FILE" "$R2_FILE" | tee "$REPORT"
