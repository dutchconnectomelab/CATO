#!/bin/bash
#
# (Example) preprocessing script for functional MRI data:
#   1. Performs slice timing correction using FSL tool slicetimer if variable
#      sliceTimingCorrection is TRUE.
#   2. Performs motion correction using FSL tool MCFLIRT.
#   3. Computes a rs-fMRI reference image by averaging all (motion corrected)
#      rs-fMRI frames (using FSL).
#   4. Computes the registration matrix between the rs-fMRI reference image and
#      the T1 image (using Freesurfer).
#   5. Registers the T1 parcellation to the reference rs-fMRI image (using
#      Freesurfer).

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
	cat << EOF
(Example) preprocessing script for functional MRI data:
    1. Performs slice timing correction using FSL tool slicetimer if
       variable sliceTimingCorrection is TRUE.
    2. Performs motion correction using FSL tool MCFLIRT.
    3. Computes a rs-fMRI reference image by averaging all (motion
       corrected) rs-fMRI frames (using FSL).
    4. Computes the registration matrix between the rs-fMRI reference
       image and the T1 image (using Freesurfer).
    5. Registers the T1 parcellation to the reference rs-fMRI image
       (using Freesurfer).
EOF
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
        --fmriFile=*)
            fmriFile=${1#*=}
            shift
        ;;
        --fmriProcessedFile=*)
            fmriProcessedFile=${1#*=}
            shift
        ;;
        --slicetimerOptions=*)
            slicetimerOptions=${1#*=}
            shift
        ;;
        --sliceTimingCorrection=*)
            sliceTimingCorrection=${1#*=}
            shift
        ;;
        --fmriInfo=*)
            if [[ ${1} == *.* ]]; then
                firstPart=${1#--fmriInfo.}
                name=${firstPart%=*}
                val=${1#*=}
                mri_convertOptions+=("-${name} ${val}")
            fi
            shift
        ;;        
        --fmriReferenceFile=*)
            fmriReferenceFile=${1#*=}
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
        --motionParametersFile=*)
            motionParametersFile=${1#*=}
            shift
        ;;     
        --freesurferDir=*)
            freesurferDir=${1#*=}
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
# Main function
#######################################
trap error ERR

parse_input "$@"

cp "$fmriFile" "$fmriProcessedFile"

# upate fmriProcessedFile header if items are missing (e.g. repetition time)
if [ ! -z "$mri_convertOptions" ]
then
    mri_convert ${mri_convertOptions[@]} "$fmriProcessedFile" "$fmriProcessedFile"
fi

tr=$(mri_info --tr  "$fmriProcessedFile" 2>/dev/null | tail -n1)
if [  $tr -eq 0 ]
then
    echo "Repetition time is missing in NifTi header (TR=0.00 msec).
    Use mri_convertOptions parameter to adjust this value to the true repetition time" >&2
    exit 1
fi

# perform slice timing correction
if [ "$sliceTimingCorrection" = true ] ; then
    slicetimer -i "$fmriProcessedFile" -o "${fmriProcessedFile/.nii.gz/_stc.nii.gz}" \
        $slicetimerOptions
    mv "${fmriProcessedFile/.nii.gz/_stc.nii.gz}" "$fmriProcessedFile"
fi

# perform motion correction
mcflirt -in "$fmriProcessedFile" -plots
mv "${fmriProcessedFile/.nii.gz/_mcf.nii.gz}" "$fmriProcessedFile"
mv "${fmriProcessedFile/.nii.gz/_mcf.par}" "$motionParametersFile"

# create reference volume
fslmaths "$fmriProcessedFile" -Tmean "$fmriReferenceFile"

# compute registration matrix
if [[ "$freesurferDir" =~ \ |\' ]] || [[ "$fmriReferenceFile" =~ \ |\' ]] || [[ "$registrationMatrixFile" =~ \ |\' ]]; then
    echo "bbregister cannot handle filenames with spaces. Variable freesurferDir, fmriReferenceFile or registrationMatrixFile contains spaces." >&2
    exit 1
fi

bbregister --s "$freesurferDir" --mov "$fmriReferenceFile" --reg "$registrationMatrixFile" \
    --bold --init-fsl
rm "$registrationMatrixFile".{mincost,param,sum,log}

# register freesurfer segmentation to reference volume
mri_label2vol --seg "${freesurferDir}/mri/aseg.mgz" --temp "$fmriReferenceFile" \
    --reg "$registrationMatrixFile" --o "$segmentationFile"

# TODO add registration of MTR weighted image

exit 0
