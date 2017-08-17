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
#                         [default: $(pwd)/ref/MM_GRCm38_p5.fasta.gz]
#   --ref-gene-map <txt>  Pass a gzip-compressed transcript-to-gene-map file
#                         [default: $(pwd)/ref/MM_GRCm38_p5_gene_map.txt.gz]
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
INPUT_REF_FASTA="$(pwd)/ref/MM_GRCm38_p5.fasta.gz"
INPUT_REF_GENE_MAP_TXT="$(pwd)/ref/MM_GRCm38_p5_gene_map.txt.gz"
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

function pull_images {
  echo '>>> pull Docker images'
  ${DC} pull
  echo
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
        pull_images && exit 0
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

echo ">>> prepare output directories
- ${OUTPUT_DIR}
  - sample
  - reference
  - summary
"
SAMPLE_DIR="${OUTPUT_DIR}/sample"
REF_DIR="${OUTPUT_DIR}/ref"
SUMMARY_DIR="${OUTPUT_DIR}/summary"
ls ${INPUT_DIR} | xargs -I {} sh -c "[[ -d ${SAMPLE_DIR}/{} ]] || mkdir -p ${SAMPLE_DIR}/{}"
[[ -d "${REF_DIR}" ]] || mkdir ${REF_DIR}
[[ -d "${SUMMARY_DIR}" ]] || mkdir ${SUMMARY_DIR}

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

# QC checks
if [[ ${NO_QC} -eq 0 ]]; then
  echo '>>> execute QC checks using FastQC'
  for s in $(ls ${SAMPLE_DIR}); do
    echo ">>>>>> sample: ${s}"
    for fq in $(ls ${INPUT_DIR}/${s}/*.fastq.gz); do
      ${DC_RUN} \
        -v ${INPUT_DIR}/${s}:/input:ro \
        -v ${SAMPLE_DIR}/${s}:/output \
        fastqc \
        --threads ${N_THREAD} \
        --nogroup \
        --outdir /output \
        /input/$(basename ${fq})
    done
  done
  echo
fi

# preparations for RSEM
RSEM_PREP_LOG="${REF_DIR}/rsem_prep_ref.log"
if [[ ! -f ${RSEM_PREP_LOG} ]]; then
  echo '>>> prepare RSEM reference files'
  REF_TAG="$(basename ${INPUT_REF_FASTA} | awk -F '.fasta' '{print $1}')"
  REF_FASTA="${REF_DIR}/${REF_TAG}.fasta"
  REF_GENE_MAP_TXT="${REF_DIR}/$(basename ${INPUT_REF_GENE_MAP_TXT} | awk -F '.txt' '{print $1".txt"}')"
  if [[ ! -f "${REF_FASTA}" ]]; then
    case "$(basename ${INPUT_REF_FASTA} | awk -F '.fasta' '{print $NF}')" in
      '.gz' )
        unpigz -p ${N_THREAD} -c ${INPUT_REF_FASTA} > ${REF_FASTA}
        ;;
      '' )
        cp ${INPUT_REF_FASTA} ${REF_FASTA}
        ;;
    esac
  fi
  if [[ ! -f "${REF_GENE_MAP_TXT}" ]]; then
    case "$(basename ${INPUT_REF_GENE_MAP_TXT} | awk -F '.txt' '{print $NF}')" in
      '.gz' )
        unpigz -p ${N_THREAD} -c ${INPUT_REF_GENE_MAP_TXT} > ${REF_GENE_MAP_TXT}
        ;;
      '' )
        cp ${INPUT_REF_GENE_MAP_TXT} ${REF_GENE_MAP_TXT}
        ;;
    esac
  fi
  ${DC_RUN} \
    -v ${REF_DIR}:/ref:ro \
    -v ${REF_RSEM_DIR}:/ref_rsem \
    --entrypoint rsem-prepare-reference \
    rsem \
    -p ${N_THREAD} \
    --transcript-to-gene-map /ref/$(basename ${REF_GENE_MAP_TXT}) \
    --bowtie2 \
    /ref/$(basename ${REF_FASTA_TXT}) \
    /ref_rsem/${REF_TAG} \
    > ${RSEM_PREP_LOG} 2>&1
  echo
fi
