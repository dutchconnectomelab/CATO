#!/bin/bash
#
# (Example) preprocessing script that:
# 1. Computes a DWI reference image based on the corrected diffusion-unweighted (b0) volumes.
# 2. Computes the registration matrix between DWI reference image and the T1 image using bbregister.
# 3. Registers the Freesurfer segmentation to the DWI reference image.

# https://google.github.io/styleguide/shellguide.html

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
# Usage
#######################################
usage()
{
    echo "(Example) preprocessing script that:"
    echo "1. Computes a DWI reference image based on the corrected diffusion-unweighted (b0) volumes."
    echo "2. Computes the registration matrix between DWI reference image and the T1 image using bbregister."
    echo "3. Registers the Freesurfer segmentation to the DWI reference image."
}

#######################################
# Parse input
#######################################
parse_input()
{
    [[ $# -eq 0 ]] && set -- "--help"
	while [ -n "$1" ]; do
    shopt -s nocasematch
    case "$1" in
   		--b0Scans=*)
        	b0Scans=${1#*=}
        	shift 		
        ;;
   		--dwiProcessedFile=*)
        	dwiProcessedFile=${1#*=}
        	shift 		
        ;;
   		--dwiFile=*)
        	dwiFile=${1#*=}
        	shift 		
        ;;
        --dwiReferenceFile=*)
        	dwiReferenceFile=${1#*=}
        	shift 		
        ;;
        --freesurferDir=*)
        	freesurferDir=${1#*=}
        	shift 		
        ;;
        --registrationMatrixFile=*)
        	registrationMatrixFile=${1#*=}
        	shift 		
        ;;        
        --segmentationFile=*)
        	segmentationFile=${1#*=}
        	shift 		
        ;;   
        -? | -h | -help | --help)
            usage && exit
        ;;             
        --*=*)
            inputParam=${1:2};
            # echo Not used: $inputParam
            shift 		
        ;;   
        *)
            error "Unkown input argument '$1'"
            shift 		
        ;;   			
    esac
done
}

#######################################
# Main function.
#######################################
trap error ERR

parse_input "$@"


# TODO: include conversion MNC to NIFTI
cp "$dwiFile" "$dwiProcessedFile"

# Load b0-scans list and parse them into an array
b0Scans=(${b0Scans//,/ })

# use b0 scan(s) to create reference volume
for i in ${b0Scans[@]}; do
	fslroi "$dwiProcessedFile" "${dwiReferenceFile/.nii.gz/_indv_$i.nii.gz}" $((i - 1)) 1
done
fslmerge -t "$dwiReferenceFile" "${dwiReferenceFile/.nii.gz/_indv_}"*.nii.gz
fslmaths "$dwiReferenceFile" -Tmean "$dwiReferenceFile"
rm "${dwiReferenceFile/.nii.gz/_indv_}"*.nii.gz

# compute registration matrix
bbregister --s "$freesurferDir" --mov "$dwiReferenceFile" --reg "$registrationMatrixFile" \
    --dti --init-fsl
rm -f "`basename $registrationMatrixFile .dat`".{dat.mincost,dat.param,dat.sum,dat.log,log}

# register freesurfer segmentation to reference volume
mri_label2vol --seg "${freesurferDir}/mri/aseg.mgz" --temp "$dwiReferenceFile" \
    --reg "$registrationMatrixFile" --o "$segmentationFile"

# TODO add registration of MTR weighted image

exit 0
