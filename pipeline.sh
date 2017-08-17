#!/usr/bin/env bash
#
# Usage:  pipeline.sh [ -h | --help | -v | --version ]
#         pipeline.sh pull
#         pipeline.sh run [--thread <int>] [--no-qc] [--input <dir>]
#                         [--output <dir>]
#
# Description:
#   Docker-based TPM calculating pipeline for single-cell RNA-seq
#
# Arguments:
#   -h, --help      Print usage
#   -v, --version   Print version information
#   --thread <int>  Limit multithreading
#   --no-qc         Skip QC checks with FastQC
#   --input <dir>   Pass an absolute path to input directory
#                   [default: $(pwd)/input]
#   --output <dir>  Pass an absolute path to output directoryr
#                   [default: $(pwd)/output]
#   pull            Pull required Docker images
#   run             Run the pipeline

set -e

[[ "${1}" = '--debug' ]] && set -x && shift 1

REPO_NAME='single-cell-tpm-calculator'
SCRIPT_NAME='pipeline.sh'
SCRIPT_VERSION='v0.0.1'
SCRIPT_ROOT="$(dirname ${0})"
SCRIPT_PATH="${SCRIPT_ROOT}/$(basename ${0})"
DC="docker-compose -f ${SCRIPT_ROOT}/docker-compose.yml"
DC_RUN="${DC} run --rm -u $(id -u):$(id -g)"
INPUT_DIR="$(pwd)/input"
OUTPUT_DIR="$(pwd)/output"
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
  -  cell
  -  reference
  -  summary
"
CELL_DIR="${OUTPUT_DIR}/cell"
REFERENCE_DIR="${OUTPUT_DIR}/reference"
SUMMARY_DIR="${OUTPUT_DIR}/summary"
ls ${INPUT_DIR} | xargs -I {} sh -c "[[ -d ${CELL_DIR}/{} ]] || mkdir -p ${CELL_DIR}/{}"
[[ -d "${REFERENCE_DIR}" ]] || mkdir ${REFERENCE_DIR}
[[ -d "${SUMMARY_DIR}" ]] || mkdir ${SUMMARY_DIR}

echo '>>> print software versions'
v_txt="${OUTPUT_DIR}/versions.txt"
print_version | tee ${v_txt}
echo -n 'fastqc: ' | tee -a ${v_txt}
${DC_RUN} fastqc --version | tee -a ${v_txt}
echo -n 'prinseq: ' | tee -a ${v_txt}
${DC_RUN} prinseq --version | tee -a ${v_txt}
echo -n 'rsem: ' | tee -a ${v_txt}
${DC_RUN} rsem --version | tee -a ${v_txt}
echo

if [[ ${NO_QC} -eq 0 ]]; then
  echo '>>> execute QC checks using FastQC'
  for sc in $(ls ${CELL_DIR}); do
    echo ">>>>>> cell: ${sc}"
    for fq in $(ls ${INPUT_DIR}/${sc}/*.fastq.gz); do
      ${DC_RUN} \
        -v ${INPUT_DIR}/${sc}:/sc_input:ro \
        -v ${CELL_DIR}/${sc}:/sc_output \
        fastqc \
        --threads ${N_THREAD} \
        --nogroup \
        --outdir /sc_output \
        /sc_input/$(basename ${fq})
    done
  done
  echo
fi
