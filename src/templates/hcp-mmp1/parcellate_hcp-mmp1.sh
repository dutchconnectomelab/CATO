#!/bin/bash
#
# Parcellate hcp-mmp1 atlas

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
        --forceFreesurferOverwrite=*)
            forceFreesurferOverwrite=${1#*=}
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

TEMPLATE=hcp-mmp1

templateDir=$(dirname "$0")
subjectDir=$(pwd)
freesurferDir=${subjectDir}/$freesurferDir
subjectFS=$(basename "${freesurferDir}") # FreeSurfer subject

# Prepare Freesurfer
SUBJECTS_DIR=$(dirname "${freesurferDir}")
cd "${SUBJECTS_DIR}"

if [ ! -f "$subjectFS/label/rh.${TEMPLATE}.annot" ] || [ "$forceFreesurferOverwrite" = true ] ; then

    # Create LH and RH annotation files
    ln -s "$FREESURFER_HOME/subjects/fsaverage" .
    mri_surf2surf --srcsubject "fsaverage" --trgsubject "${subjectFS}" --hemi lh   \
        --sval-annot "${templateDir}/lh.${TEMPLATE}.annot"  --tval "$subjectFS/label/lh.${TEMPLATE}.annot"
    mri_surf2surf --srcsubject "fsaverage" --trgsubject "${subjectFS}" --hemi rh   \
        --sval-annot "${templateDir}/rh.${TEMPLATE}.annot"  --tval "$subjectFS/label/rh.${TEMPLATE}.annot"
    unlink "fsaverage"

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
