#!/bin/bash

#########################
# Clustering algorithms #
#########################

# This script crates a PBS job for our snakemake clustering pipeline

# set -euo pipefail

function usage () {
    cat >&2 <<EOF

USAGE: ${0} [-y] [options]
  -y <config file> : Path to the YAML config file. Required.
  -s Submit job.
  -d Job ID after which, if successfully finishes, it will start running.
  -v Verbose.
  -h Print the usage info.

EOF
}
# initial : makes this loop silent and now requires '?)'
# ${opt} is each option and ${OPTARG} its the argumet (if a colon is there ${opt}:)
SUBMIT=FALSE
VERBOSE=FALSE
while getopts ":y:sd:vh" opt; do
  case ${opt} in
    y) CONFIG_FILE=${OPTARG};;
    s) SUBMIT=TRUE;;
    d) DEPEND=${OPTARG:-""};;
    v) VERBOSE=TRUE;;
    h) usage; exit 1;;
    \?) echo "No -${OPTARG} argument found."; usage; exit 1;;
  esac
done
if [[ ${OPTIND} -eq 1 ]] ; then
    usage; exit 1
fi

#### Parameters #### -----------------------------------------------------------
function read_yaml(){
  sed 's/#.*//g' ${1} | grep ${2}: | sed 's/.*:[^:\/\/]//; s/\"//g'
}
OUTPUT_DIR="$(read_yaml ${CONFIG_FILE} output_dir)"
PROJECT_NAME="$(read_yaml ${CONFIG_FILE} project_name)"
TOOLS="$(read_yaml ${CONFIG_FILE} tool | sed -E 's/\[|\]//g; s/,/_/g; s/ {1,}//g')"
OUTPUT_DIR="${OUTPUT_DIR%/}/${PROJECT_NAME}"

if grep -q 'pipeline:' "${CONFIG_FILE}"; then
  PIPELINE_DIR=$(grep 'pipeline:' "${CONFIG_FILE}" | awk '{print $2}' | sed 's/\"//g')
else
  PIPELINE_DIR=$(dirname "${0}")
fi
if [[ ! -s "${PIPELINE_DIR}/cluster.json" ]]; then PIPELINE_DIR=$(dirname "${0}"); fi
PIPELINE_DIR="${PIPELINE_DIR%/}"
CLUST_ENVIRON="$(read_yaml ${CONFIG_FILE} environment)"
CLUST_CONFIG="$(read_yaml ${CONFIG_FILE} cluster_config)"
if [[ ! -s "${CLUST_CONFIG}" ]]; then
  CLUST_CONFIG="${PIPELINE_DIR}/cluster.json"
fi

echo ' '
echo -e "\033[0;36m**** Vijay Lab - LJI 2020\033[0m"

echo -e "\033[0;36m------------------------------- PRESENTING PARAMETERS -------------------------------\033[0m"
echo "Configuration file: ${CONFIG_FILE}"
echo "Output path: ${OUTPUT_DIR}"
echo "Pipeline: ${PIPELINE_DIR}"
echo "Cluster configuration: ${CLUST_CONFIG}"
echo -e "\033[0;36m------------------------------- --------------------- -------------------------------\033[0m"

if [[ ! -d "${OUTPUT_DIR}" ]]; then mkdir --parents "${OUTPUT_DIR}"; fi
if [[ ! -d "${OUTPUT_DIR}/scripts" ]]; then mkdir "${OUTPUT_DIR}/scripts"; fi

cp ${CONFIG_FILE} ${OUTPUT_DIR}/config.yaml # migh consider running a check

JOBFILE="${OUTPUT_DIR}/scripts/clump_${PROJECT_NAME}_${TOOLS}"
echo "Job file: ${JOBFILE}.sh"
rm ${OUTPUT_DIR}/scripts/clump_${PROJECT_NAME}*txt 2> /dev/null

wget -q https://raw.githubusercontent.com/vijaybioinfo/cellranger_wrappeR/main/routine_template.sh -O ${JOBFILE}.sh

sed -i 's|{cellranger}|clustering|' ${JOBFILE}.sh
sed -i 's|{username}|'"${USER}"'|g' ${JOBFILE}.sh
sed -i 's|{sampleid}|'"${PROJECT_NAME}"'|g' ${JOBFILE}.sh
sed -i 's|\/\.\.||g' ${JOBFILE}.sh
sed -i 's|{routine_pbs}|clump|' ${JOBFILE}.sh
sed -i 's|{outpath}|'"${OUTPUT_DIR}"'|g' ${JOBFILE}.sh
sed -i 's|cp ${PROJ.*|cp -r ${PROJDIR}/. ./|g' ${JOBFILE}.sh # to copy everything to scratch
sed -i 's|cp -R ./.*${PROJ.*|cp -r . ${PROJDIR}/|g' ${JOBFILE}.sh # copy from scratch
echo "Pushing critical lines..."
if [[ "${CLUST_ENVIRON}" != "" ]] && [[ "$(which conda | wc -l)" == "1" ]]; then
  echo "Environment: $(conda env list | grep cluster | awk '{print $2 $3}')"
  sed -i 's|# {after_copy}|source activate '"${CLUST_ENVIRON}"'; conda env list|g' ${JOBFILE}.sh
  sed -i 's|# {pre_routine}|snakemake --snakefile '"${PIPELINE_DIR}"'/Snakefile --dag \| dot -Tpdf > '"${OUTPUT_DIR}"'/_results_outline.pdf|g' ${JOBFILE}.sh
fi
sed -i 's|{routine_params}|snakemake --jobs 100 --latency-wait 100 --cluster-config '"${CLUST_CONFIG}"' --snakefile '"${PIPELINE_DIR}"'/Snakefile --cluster "qsub -l {cluster.walltime} -l {cluster.cores} -l {cluster.memory} -m n -q default -e '"${OUTPUT_DIR}"'/scripts/ -o '"${OUTPUT_DIR}"'/scripts/" --jobname "clump.{rulename}.{jobid}" --stats '"${OUTPUT_DIR}"'/scripts/snakemake.stats >\& '"${OUTPUT_DIR}"'/scripts/snakemake.log|' ${JOBFILE}.sh

sed -i 's|{walltime}|36:00:00|g' ${JOBFILE}.sh
sed -i 's|{nodes}|1|g' ${JOBFILE}.sh
sed -i 's|{ppn}|1|g' ${JOBFILE}.sh
sed -i 's|{mem}|8gb|g' ${JOBFILE}.sh

if echo "${SUBMIT}" | grep -qE "TRUE|^yes$|^y$"; then
  echo "Check it out"; exit
fi
if [[ "${DEPEND}" != "" ]]; then
  DEPEND=$(qsub -W depend=afterok:${DEPEND} ${JOBFILE}.sh)
else
  CID=$(qsub ${JOBFILE}.sh)
fi; CID=$(echo "${CID}" | sed 's/\..*//')
echo "Job ID: ${CID}";
echo
