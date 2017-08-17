#!/usr/bin/env bash
#
# Usage:  pipeline.sh [ -h | --help | -v | --version ]
#         pipeline.sh pull
#         pipeline.sh run [--thread <int>] [--no-qc] [--ref <fasta>]
#                         [--ref-gene-map <txt>] [--input <dir>]
#                         [--output <dir>]
#
# Description:
#   Docker-based TPM calculating pipeline for RNA-seq
#
# Arguments:
#   -h, --help            Print usage
#   -v, --version         Print version information
#   --thread <int>        Limit multithreading
#   --no-qc               Skip QC checks with FastQC
#   --ref <fasta>         Pass a gzip-compressed reference FASTA file
#                         [default: $(pwd)/ref/GRCm38_p5.fasta.gz]
#   --ref-gene-map <txt>  Pass a gzip-compressed transcript-to-gene-map file
#                         [default: $(pwd)/ref/GRCm38_p5_gene_map.txt.gz]
#   --input <dir>         Pass an absolute path to input directory
#                         [default: $(pwd)/input]
#   --output <dir>        Pass an absolute path to output directoryr
#                         [default: $(pwd)/output]
#   pull                  Pull required Docker images
#   run                   Run the pipeline

set -e

[[ "${1}" = '--debug' ]] && set -x && shift 1

REPO_NAME='rna-tpm-calculator'
SCRIPT_NAME='pipeline.sh'
SCRIPT_VERSION='v0.0.1'
SCRIPT_ROOT="$(dirname ${0})"
SCRIPT_PATH="${SCRIPT_ROOT}/$(basename ${0})"
DC="docker-compose -f ${SCRIPT_ROOT}/docker-compose.yml"
DC_RUN="${DC} run --rm -u $(id -u):$(id -g)"
INPUT_DIR="$(pwd)/input"
OUTPUT_DIR="$(pwd)/output"
INPUT_REF_FASTA="$(pwd)/ref/GRCm38_p5.fasta.gz"
INPUT_REF_GENE_MAP_TXT="$(pwd)/ref/GRCm38_p5_gene_map.txt.gz"
RUN_FLAG=0
NO_QC=0

function print_version {
  echo "${REPO_NAME}: ${SCRIPT_NAME} ${SCRIPT_VERSION}"
}

function print_usage {
  sed -ne '1,2d; /^#/!q; s/^#$/# /; s/^# //p;' ${SCRIPT_PATH}
}

function abort {
  echo " ${SCRIPT_NAME}: ${*}" >&2
  exit 1
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
        INPUT_REF_FASTA="$(dirname ${2})/$(basename ${2})" && shift 2
        ;;
      '--ref-gene-map' )
        INPUT_REF_GENE_MAP_TXT="$(dirname ${2})/$(basename ${2})" && shift 2
        ;;
      '--input' )
        INPUT_DIR="$(dirname ${2})/$(basename ${2})" && shift 2
        ;;
      '--output' )
        OUTPUT_DIR="$(dirname ${2})/$(basename ${2})" && shift 2
        ;;
      'pull' )
        echo '>>> pull images' && ${DC} pull && exit 0
        ;;
      'run' )
        RUN_FLAG=1 && shift 1
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
PIGZ="pigz -p ${N_THREAD}"
UNPIGZ="unpigz -p ${N_THREAD}"

### 1.  Preparation
echo ">>> prepare output directories
- ${OUTPUT_DIR}
  - sample
  - ref
  - summary"
SAMPLE_DIR="${OUTPUT_DIR}/sample"
REF_DIR="${OUTPUT_DIR}/ref"
SUMMARY_DIR="${OUTPUT_DIR}/summary"
ls ${INPUT_DIR} | xargs -I {} bash -c "\
  [[ -d ${SAMPLE_DIR}/{} ]] || mkdir -p ${SAMPLE_DIR}/{}
"
[[ -d "${REF_DIR}" ]] || mkdir ${REF_DIR}
[[ -d "${SUMMARY_DIR}" ]] || mkdir ${SUMMARY_DIR}
echo

echo '>>> print software versions'
VERSIONS_TXT="${OUTPUT_DIR}/versions.txt"
print_version | tee ${VERSIONS_TXT}
echo -n 'fastqc: ' | tee -a ${VERSIONS_TXT}
${DC_RUN} fastqc --version | tee -a ${VERSIONS_TXT}
echo -n 'prinseq: ' | tee -a ${VERSIONS_TXT}
${DC_RUN} prinseq --version | tee -a ${VERSIONS_TXT}
echo -n 'rsem: ' | tee -a ${VERSIONS_TXT}
${DC_RUN} rsem --version | tee -a ${VERSIONS_TXT}
echo

echo '>>> concatenate fastq files'
ls ${SAMPLE_DIR} | xargs -P ${N_THREAD} -I {} bash -c "\
  cat ${INPUT_DIR}/{}/*.fastq.gz > ${SAMPLE_DIR}/{}/raw.fastq.gz
"
echo

REF_TAG="$(basename ${INPUT_REF_FASTA} | awk -F '.fasta' '{print $1}')"
RSEM_PREP_LOG="${REF_DIR}/rsem_prepare_reference.log"
if [[ ! -f ${RSEM_PREP_LOG} ]]; then
  echo '>>> prepare RSEM reference files'
  mkdir "${REF_DIR}/rsem"
  ${UNPIGZ} -c ${INPUT_REF_FASTA} > ${REF_DIR}/${REF_TAG}.fasta
  ${UNPIGZ} -c ${INPUT_REF_GENE_MAP_TXT} > ${REF_DIR}/${REF_TAG}_gene_map.txt
  ${DC_RUN} -v ${REF_DIR}:/ref:ro -v ${REF_DIR}/rsem:/ref_rsem \
    --entrypoint rsem-prepare-reference \
    rsem \
    -p ${N_THREAD} \
    --bowtie2 \
    --transcript-to-gene-map /ref/${REF_TAG}_gene_map.txt \
    /ref/${REF_TAG}.fasta \
    /ref_rsem/${REF_TAG} \
    > ${RSEM_PREP_LOG} 2>&1
  ${PIGZ} ${REF_DIR}/${REF_TAG}.fasta ${REF_DIR}/${REF_TAG}_gene_map.txt
  echo
fi


### 2.  QC checks
if [[ ${NO_QC} -eq 0 ]]; then
  echo '>>> execute QC checks using FastQC'
  ls ${SAMPLE_DIR} | xargs -I {} bash -c "\
    ${DC_RUN} -v ${SAMPLE_DIR}/{}:/output \
    fastqc \
    --threads ${N_THREAD} \
    --nogroup \
    --outdir /output \
    /output/raw.fastq.gz \
    > ${SAMPLE_DIR}/{}/fastqc.log 2>&1
  "
  echo
fi


### 3.  Trimming and filtering
echo '>>> trim and filter sequences'
find ${SAMPLE_DIR} -maxdepth 2 -type f -name 'raw.fastq.gz' | xargs ${UNPIGZ}
ls ${SAMPLE_DIR} | xargs -P ${N_THREAD} -I {} bash -c "\
  ${DC_RUN} -v ${SAMPLE_DIR}/{}:/output \
  prinseq \
  -min_len 30 \
  -trim_tail_right 5 \
  -trim_tail_left 5 \
  -min_qual_mean 20 \
  -trim_qual_right 20 \
  -trim_qual_left 20 \
  -fastq output/raw.fastq \
  -out_good output/prinseq_good \
  -out_bad output/prinseq_bad \
  > ${SAMPLE_DIR}/{}/prinseq_lite.log 2>&1
"
find ${SAMPLE_DIR} -maxdepth 2 -type f -name '*.fastq' | xargs ${PIGZ}
echo


### 4.  Mapping
echo '>>> map reads and calculate TPM'
ls ${SAMPLE_DIR} | xargs -I {} bash -c "\
  ${DC_RUN} -v ${SAMPLE_DIR}/{}:/output -v ${REF_DIR}/rsem:/ref_rsem \
  --entrypoint rsem-calculate-expression \
  rsem \
  --bowtie2 \
  --estimate-rspd \
  -p ${N_THREAD} \
  --calc-ci \
  /output/prinseq_good.fastq.gz \
  /ref_rsem/${REF_TAG} \
  /output/rsem_${REF_TAG}_{} \
  > ${SAMPLE_DIR}/{}/rsem_calculate_expression.log 2>&1
"
echo
