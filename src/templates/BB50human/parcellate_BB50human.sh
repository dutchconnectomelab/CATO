#!/bin/bash
#
# Parcellate BB50human atlas

# set -x
set -e

#######################################
# Error handling
#######################################
error() {
    if [ -n "$1" ]; then
    echo -e "error: $1\n" >&2
    else
        read line file <<<$(caller)
        echo "An error occurred in line $line of file $file:" >&2
        sed "${line}q;d" "$file" >&2
    fi
    exit 1
}

#######################################
# Parse input
#######################################
parse_input()
{

    while [ -n "$1" ]; do
    shopt -s nocasematch
    case "$1" in
        --freesurferDir=*)
            freesurferDir=${1#*=}
            shift       
        ;;
        --referenceFile=*)
            referenceFile=${1#*=}
            shift       
        ;;
        --registrationMatrixFile=*)
            registrationMatrixFile=${1#*=}
            shift       
        ;;
        --parcellationFile=*)
            parcellationFile=${1#*=}
            shift       
        ;;            
        --forceFreesurferOVerwrite=*)
            forceFreesurferOVerwrite=${1#*=}
            shift       
        ;;
    --*=*)
        inputParam=${1:2};
        shift       
        ;;   
    *)
        error Unkown input argument: $inputParam
        shift       
        ;;              
    esac
done
}

#######################################
# Main function
#######################################
trap error ERR
parse_input "$@"

TEMPLATE=BB50human

templateDir=$(dirname "$0") # describing the directory that includes the lh.colortable.txt, rh.colortable.txt, lh.economo.gcs, rh.economo.gcs datafiles
subjectDir=$(pwd)
freesurferDir=${subjectDir}/$freesurferDir
subjectFS=$(basename "${freesurferDir}") # FreeSurfer subject

# Prepare Freesurfer
SUBJECTS_DIR=$(dirname "${freesurferDir}")
cd "${SUBJECTS_DIR}"

if [ ! -f "$subjectFS/label/rh.${TEMPLATE}.annot" ] || [ "$forceFreesurferOverwrite" = true ] ; then

	# Create LH and RH annotation files
	mris_ca_label -t "${templateDir}/${TEMPLATE}.annot.ctab" "${subjectFS}" lh "${subjectFS}/surf/lh.sphere.reg" "${templateDir}/lh.${TEMPLATE}.gcs" "${subjectFS}/label/lh.${TEMPLATE}.annot"
	mris_ca_label -t "${templateDir}/${TEMPLATE}.annot.ctab" "${subjectFS}" rh "${subjectFS}/surf/rh.sphere.reg" "${templateDir}/rh.${TEMPLATE}.gcs" "${subjectFS}/label/rh.${TEMPLATE}.annot"

	# Add cortical labels to the automatic segmentation volume (aseg)
	mri_aparc2aseg --s "${subjectFS}" --annot "${TEMPLATE}"

fi

if [ ! -f "$subjectFS/stats/rh.${TEMPLATE}.stats" ] || [ "$forceFreesurferOverwrite" = true ] ; then

	# Create anatomical stat files
	mris_anatomical_stats -a "${subjectFS}/label/lh.${TEMPLATE}.annot" -f "${subjectFS}/stats/lh.${TEMPLATE}.stats" "${subjectFS}" lh
	mris_anatomical_stats -a "${subjectFS}/label/rh.${TEMPLATE}.annot" -f "${subjectFS}/stats/rh.${TEMPLATE}.stats" "${subjectFS}" rh

fi


# Derive segmentation of dwi B0-reference volume
cd "${subjectDir}"
mri_label2vol --seg "${freesurferDir}/mri/${TEMPLATE}+aseg.mgz" \
              --temp "${referenceFile}" --reg "${registrationMatrixFile}" --o "${parcellationFile}"

exit 0
