#!/usr/bin/env bash
#
# Usage:  workflow.sh [ -h | --help | -v | --version ]
#         workflow.sh run [--thread <int>] [--no-qc] [--ref <fasta>]
#                         [--ref-gene-map <txt>] [--input <dir>]
#                         [--output <dir>]
#
# Description:
#   TPM calculating workflow for RNA-seq
#
# Arguments:
#   -h, --help              Print usage
#   -v, --version           Print version information
#   --thread <int>          Limit multithreading
#   --no-qc                 Skip QC checks with FastQC
#   --ref <fasta>           Pass a gzip-compressed reference FASTA file
#                           [default: ./ref/GRCm38_p5.fasta.gz]
#   --ref-gene-map <txt>    Pass a gzip-compressed transcript-to-gene-map file
#                           [default: ./ref/GRCm38_p5_gene_map.txt.gz]
#   --input <dir>           Pass an path to input directory
#                           [default: ./input]
#   --output <dir>          Pass an path to output directoryr
#                           [default: ./output]
#   run                     Run the workflow

set -e

[[ "${1}" = '--debug' ]] && set -x && shift 1

REPO_NAME='rna-tpm-calculator'
REPO_VERSION='v0.0.1'
REPO_ROOT="$(dirname ${0})"
WORKFLOW_SH="${REPO_ROOT}/$(basename ${0})"
INPUT_DIR="$(pwd)/input"
OUTPUT_DIR="$(pwd)/output"
INPUT_REF_FASTA="$(pwd)/ref/GRCm38_p5.fasta.gz"
INPUT_REF_GENE_MAP_TXT="$(pwd)/ref/GRCm38_p5_gene_map.txt.gz"
NO_QC=0

function print_version {
  echo "${REPO_NAME}: ${REPO_VERSION}"
}

function print_usage {
  sed -ne '1,2d; /^#/!q; s/^#$/# /; s/^# //p;' ${WORKFLOW_SH}
}

function abort {
  echo " ${REPO_NAME}: ${*}" >&2
  exit 1
}

function print_abspath {
  python -c "import os; print(os.path.abspath('${*}'))"
}

function count_cpu_cores {
  case "${OSTYPE}" in
    darwin*)
      sysctl -n hw.ncpu
      ;;
    linux*)
      grep -ce '^processor' /proc/cpuinfo
      ;;
  esac
}

N_THREAD=$(count_cpu_cores)

if [[ -n "${1}" ]]; then
  while [[ -n "${1}" ]]; do
    case "${1}" in
      '-v' | '--version' )
        print_version && exit 0
        ;;
      '-h' | '--help' )
        print_usage && exit 0
        ;;
      '--thread' )
        N_THREAD=${2} && shift 2
        ;;
      '--no-qc' )
        NO_QC=1 && shift 1
        ;;
      '--ref' )
        INPUT_REF_FASTA="$(print_abspath ${2})" && shift 2
        ;;
      '--ref-gene-map' )
        INPUT_REF_GENE_MAP_TXT="$(print_abspath ${2})" && shift 2
        ;;
      '--input' )
        INPUT_DIR="$(print_abspath ${2})" && shift 2
        ;;
      '--output' )
        OUTPUT_DIR="$(print_abspath ${2})" && shift 2
        ;;
      'run' )
        shift 1
        ;;
      * )
        abort "invalid argument \`${1}\`"
        ;;
    esac
  done
else
  print_usage && exit 0
fi

set -u

PGZ="pigz -p ${N_THREAD}"
CREATE_MATRIX_R="$(print_abspath ${REPO_ROOT}/create_matrix.R)"
SAMPLE_DIR="${OUTPUT_DIR}/sample"
REF_DIR="${OUTPUT_DIR}/ref"
REF_RSEM_DIR="${REF_DIR}/rsem"
SUMMARY_DIR="${OUTPUT_DIR}/summary"
READ_COUNT_CSV="${SUMMARY_DIR}/read_count.csv"
COMPLETED_LOG="${OUTPUT_DIR}/completed.log"
VER_TXT="${OUTPUT_DIR}/versions.txt"
REF_TAG="$(basename ${INPUT_REF_FASTA} | awk -F '.fasta' '{print $1}')"
WORKFLOW=(
  '>>> make output directories'
  '>>> concatenate fastq files'
  '>>> prepare RSEM reference files'
  '>>> execute QC checks using FastQC'
  '>>> trim and filter sequences'
  '>>> map reads and calculate TPM'
)

function is_not_completed {
  [[ $(grep -ce "^completed - .* - ${*}$" ${COMPLETED_LOG}) -eq 0 ]]
}

function echo_completed {
  echo "completed - $(date) - ${*}" | tee -a ${COMPLETED_LOG}
}


# 0.  make output directories
if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "${WORKFLOW[0]}
  - ${OUTPUT_DIR}
    - sample
    - ref
    - summary"
  mkdir -p ${SAMPLE_DIR} ${REF_RSEM_DIR} ${SUMMARY_DIR}
  echo "[${REPO_NAME}]" | tee ${VER_TXT} && print_version | tee -a ${VER_TXT}
  echo_completed ${WORKFLOW[0]}
elif [[ ! -f "${COMPLETED_LOG}" ]]; then
  touch ${COMPLETED_LOG}
fi
echo 'id,src,type,count' > ${READ_COUNT_CSV}


# 1.  concatenate fastq files
if $(is_not_completed ${WORKFLOW[1]}); then
  echo "${WORKFLOW[1]}"
  ls -L ${INPUT_DIR} | xargs -P ${N_THREAD} -I {} bash -c \
    "mkdir ${SAMPLE_DIR}/{} && cat ${INPUT_DIR}/{}/*.fastq.gz > ${SAMPLE_DIR}/{}/raw.fastq.gz"
  find ${SAMPLE_DIR} -name 'raw.fastq.gz' \
    | xargs -P ${N_THREAD} zgrep -ce '^@[A-Z0-9\-]\+:' \
    | sed -e 's/^.*\/\([^\/]\+\)\/raw\.fastq\.gz:\([0-9]\+\)$/\1,fastq,raw,\2/' \
    | tee -a ${READ_COUNT_CSV}
  echo_completed ${WORKFLOW[1]}
fi


# 2.  prepare RSEM reference files
if $(is_not_completed ${WORKFLOW[2]}); then
  echo "${WORKFLOW[2]}"
  echo '[bowtie2]' | tee -a ${VER_TXT} && bowtie2 --version | tee -a ${VER_TXT}
  echo '[rsem]' | tee -a ${VER_TXT} && rsem-calculate-expression --version | tee -a ${VER_TXT}
  ${PGZ} -dc ${INPUT_REF_FASTA} > ${REF_DIR}/${REF_TAG}.fasta
  ${PGZ} -dc ${INPUT_REF_GENE_MAP_TXT} > ${REF_DIR}/${REF_TAG}_gene_map.txt
  rsem-prepare-reference \
    -p ${N_THREAD} \
    --bowtie2 \
    --transcript-to-gene-map ${REF_DIR}/${REF_TAG}_gene_map.txt \
    ${REF_DIR}/${REF_TAG}.fasta \
    ${REF_RSEM_DIR}/${REF_TAG} \
    > ${REF_RSEM_DIR}/rsem_prepare_reference.log 2>&1
  ${PGZ} ${REF_DIR}/${REF_TAG}.fasta ${REF_DIR}/${REF_TAG}_gene_map.txt
  echo_completed ${WORKFLOW[2]}
fi


# 3.  execute QC checks using FastQC
if [[ ${NO_QC} -ne 0 ]] && $(is_not_completed ${WORKFLOW[3]}); then
  echo "${WORKFLOW[3]}"
  echo '[fastqc]' | tee -a ${VER_TXT} && fastqc --version | tee -a ${VER_TXT}
  for s in $(ls -L ${SAMPLE_DIR}); do
    fastqc \
      --threads ${N_THREAD} \
      --nogroup \
      --outdir ${SAMPLE_DIR}/${s} \
      ${SAMPLE_DIR}/${s}/raw.fastq.gz \
      > ${SAMPLE_DIR}/${s}/fastqc.log 2>&1
  done
  echo_completed ${WORKFLOW[3]}
fi


# 4.  trim and filter sequences
if $(is_not_completed ${WORKFLOW[4]}); then
  echo "${WORKFLOW[4]}"
  echo '[prinseq]' | tee -a ${VER_TXT} && prinseq-lite.pl --version | tee -a ${VER_TXT}
  ls -L ${SAMPLE_DIR} | xargs -P ${N_THREAD} -I {} bash -c "\
    gzip -dc ${SAMPLE_DIR}/{}/raw.fastq.gz | perl /usr/local/src/prinseq/prinseq-lite.pl \
    -min_len 30 -trim_tail_right 5 -trim_tail_left 5 \
    -min_qual_mean 20 -trim_qual_right 20 -trim_qual_left 20 \
    -fastq stdin \
    -out_good ${SAMPLE_DIR}/{}/prinseq_good \
    -out_bad ${SAMPLE_DIR}/{}/prinseq_bad \
    > ${SAMPLE_DIR}/{}/prinseq_lite.log 2>&1"
  find ${SAMPLE_DIR} -name 'prinseq_bad.fastq' | xargs ${PGZ}
  find ${SAMPLE_DIR} -name 'prinseq_good.fastq' \
    | xargs -P ${N_THREAD} zgrep -ce '^@[A-Z0-9\-]\+:' \
    | sed -e 's/^.*\/\([^\/]\+\)\/prinseq_good\.fastq:\([0-9]\+\)$/\1,fastq,qc,\2/' \
    | tee -a ${READ_COUNT_CSV}
  echo_completed ${WORKFLOW[4]}
fi


# 5.  map reads and calculate TPM
if $(is_not_completed ${WORKFLOW[5]}); then
  echo "${WORKFLOW[5]}"
  for s in $(ls -L ${SAMPLE_DIR}); do
    rsem-calculate-expression \
      --bowtie2 \
      --estimate-rspd \
      -p ${N_THREAD} \
      --calc-ci \
      ${SAMPLE_DIR}/${s}/prinseq_good.fastq \
      ${REF_RSEM_DIR}/${REF_TAG} \
      ${SAMPLE_DIR}/${s}/rsem_aligned \
      > ${SAMPLE_DIR}/${s}/rsem_calculate_expression.log 2>&1
  done
  find ${SAMPLE_DIR} -name 'prinseq_good.fastq' | xargs ${PGZ}
  find ${SAMPLE_DIR} -name 'rsem_aligned.transcript.bam' \
    | sed -e 's/\/rsem_aligned.transcript.bam$//' \
    | xargs -P ${N_THREAD} -I {} bash -c \
    'samtools flagstat {}/rsem_aligned.transcript.bam > {}/samtools_flagstat.txt'
  find ${SAMPLE_DIR} -name 'samtools_flagstat.txt' \
    | sed -e 's/^.*\/\([^\/]\+\)\/samtools_flagstat\.txt$/\1/' \
    | xargs -P ${N_THREAD} -I {} awk '{
      if (NR == 1) {
        print "{},bam,total,"$1
      } else if (NR == 2) {
        print "{},bam,secondary,"$1
      } else if (NR == 5) {
        print "{},bam,mapped,"$1
      }
    }' ${SAMPLE_DIR}/{}/samtools_flagstat.txt \
    | tee -a ${READ_COUNT_CSV}
  Rscript ${CREATE_MATRIX_R} --output "${OUTPUT_DIR}" > ${SUMMARY_DIR}/create_matrix.log
  echo_completed ${WORKFLOW[5]}
fi
