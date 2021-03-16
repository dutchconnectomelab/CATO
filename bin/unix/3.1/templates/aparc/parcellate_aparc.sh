#!/bin/bash
#
# Parcellate aparc atlas

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

subjectDir=$(pwd)
freesurferDir=${subjectDir}/${freesurferDir}

# Create parcellation
cd "${subjectDir}"
mri_label2vol --seg "${freesurferDir}/mri/aparc+aseg.mgz" \
              --temp "${referenceFile}" --reg "${registrationMatrixFile}" --o "${parcellationFile}"


exit 0
