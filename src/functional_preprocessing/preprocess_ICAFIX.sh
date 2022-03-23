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
#   6. Perform FSL tool melodic.
#   7. Perform FSL tool FIX.

# Needs configuration file parameters, e.g.:
#   "highPass":2000, 
#   "fixDir":"/usr/local/fix", 
#   "trainingData":"/usr/local/fix/training_files/standard.RData", 
#   "fixThreshold":5, 
#   "doMotionRegression":true, 
#   "T1File":"T1/SUBJECT_T1.nii.gz"    

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
    6. Perform FSL tool melodic.
    7. Perform FSL tool FIX.

Needs configuration file parameters, e.g.:
      "highPass":2000, 
      "fixDir":"/usr/local/fix", 
      "trainingData":"/usr/local/fix/training_files/standard.RData", 
      "fixThreshold":5, 
      "doMotionRegression":true, 
      "T1File":"T1/SUBJECT_T1.nii.gz"    
       
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
        --mri_convertOptions*)
            firstPart=${1#--mri_convertOptions.}
            name=${firstPart%=*}
            val=${1#*=}
            mri_convertOptions+=("-${name} ${val}")
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
        --highPass=*)
            highPass=${1#*=}
            shift
        ;; 
        --trainingData=*)
            trainingData=${1#*=}
            shift
        ;; 
        --fixThreshold=*)
            fixThreshold=${1#*=}
            shift
        ;; 
        --doMotionRegression=*)
            doMotionRegression=${1#*=}
            shift
        ;; 
        --fixDir=*)
            fixDir=${1#*=}
            shift
        ;;        
        --T1File=*)
            T1File=${1#*=}
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

# Upate fmriProcessedFile header if items are missing (e.g. RT)
if [ ! -z "$mri_convertOptions" ]
then
    mri_convert ${mri_convertOptions[@]} "$fmriProcessedFile" "$fmriProcessedFile"
fi

tr=$(mri_info --tr  $fmriProcessedFile 2>/dev/null | tail -n1)
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

# create reference volume
fslmaths "$fmriProcessedFile" -Tmean "$fmriReferenceFile"

# compute registration matrix
if [[ "$freesurferDir" =~ \ |\' ]] \
    || [[ "$fmriReferenceFile" =~ \ |\' ]] \
    || [[ "$registrationMatrixFile" =~ \ |\' ]]; then
    echo "bbregister cannot handle filenames with spaces.
    Variable freesurferDir, fmriReferenceFile or registrationMatrixFile contains spaces." >&2
    exit 1
fi

bbregister --s "$freesurferDir" --mov "$fmriReferenceFile" --reg "$registrationMatrixFile" \
    --bold --init-fsl
rm "$registrationMatrixFile".{mincost,param,sum,log}

# register freesurfer segmentation to reference volume
mri_label2vol --seg "${freesurferDir}/mri/aseg.mgz" --temp "$fmriReferenceFile" \
    --reg "$registrationMatrixFile" --o "$segmentationFile"

# perform motion correction
echo "Perform motion correction"
mcflirt -in "$fmriProcessedFile" -plots
mv "${fmriProcessedFile/.nii.gz/_mcf.nii.gz}" "$fmriProcessedFile"
mv "${fmriProcessedFile/.nii.gz/_mcf.par}" "$motionParametersFile"

## MELODIC
# High-pass threshold
fmrihp=${fmriProcessedFile}_hp${highPass}
hptr=$(echo "scale = 10; $highPass / (2 * $tr)" | bc -l)
fslmaths ${fmriProcessedFile} -Tmean ${fmrihp}
fslmaths ${fmriProcessedFile} -sub ${fmrihp} -bptf ${hptr} -1 -add ${fmrihp} ${fmrihp}

# run melodic
echo "Run melodic"
fmrihpDir="${fmrihp}_ICA_FIX"
mkdir -p ${fmrihpDir}
melodic -i ${fmrihp} -o ${fmrihpDir}/filtered_func_data.ica --nobet --report --Oall --tr=${tr}

## FIX
# Get Movement_Regressors.txt into the format expected by functionmotionconfounds.m
mkdir -p ${fmrihpDir}/mc
cat $motionParametersFile | awk '{ print $4 " " $5 " " $6 " " $1 " " $2 " " $3}' \
> ${fmrihpDir}/mc/prefiltered_func_data_mcf.par

ln -sf $(pwd)/$fmriProcessedFile ${fmrihpDir}/filtered_func_data.nii.gz
ln -sf $(pwd)/$fmriReferenceFile ${fmrihpDir}/mean_func.nii.gz  
ln -sf $(pwd)/${fmrihpDir}/filtered_func_data.ica/mask.nii.gz ${fmrihpDir}/mask.nii.gz

# make reg directory
mkdir -p ${fmrihpDir}/reg
ln -sf $(pwd)/$fmriReferenceFile ${fmrihpDir}/reg/example_func.nii.gz
ln -sf $(pwd)/$T1File ${fmrihpDir}/reg/highres.nii.gz
flirt -in $(pwd)/$T1File \
      -ref ${fmrihpDir}/reg/example_func.nii.gz \
      -omat ${fmrihpDir}/reg/highres2example_func.mat

# run fix
echo "Run fix"
if [[ ${doMotionRegression} = "true" ]]; then
    ${fixDir}/fix "$(pwd)/${fmrihpDir}" "${trainingData}" "${fixThreshold}" -m -h "${highPass}"
else
    ${fixDir}/fix "$(pwd)/${fmrihpDir}" "${trainingData}" "${fixThreshold}" -h "${highPass}"
fi

cp "$(pwd)/${fmrihpDir}/filtered_func_data_clean.nii.gz" "${fmriProcessedFile}"

exit 0
