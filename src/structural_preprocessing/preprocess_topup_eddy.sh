#!/bin/bash
#
# (Example) preprocessing script that:
# 1. Runs FSL eddy and FSL topup to correct for susceptibility induced distortions, eddy current distortions and motion artifacts in the DWI data.
# 2. Updates the b-vectors to adjust for the DWI corrections.
# 3. Computes a DWI reference image based on the corrected diffusion-unweighted (b0) volumes.
# 4. Computes the registration matrix between DWI reference image and the T1 image using bbregister.
# 5. Registers the Freesurfer segmentation to the DWI reference image.

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
    echo "1. Runs FSL eddy and FSL topup to correct for susceptibility induced distortions, eddy current distortions and motion artifacts in the DWI data."
    echo "2. Updates the b-vectors to adjust for the DWI corrections."
    echo "3. Computes a DWI reference image based on the corrected diffusion-unweighted (b0) volumes."
    echo "4. Computes the registration matrix between DWI reference image and the T1 image using bbregister."
    echo "5. Registers the Freesurfer segmentation to the DWI reference image."
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
   		--dwiFiles=*)
        	dwiFiles=${1#*=}
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
        --processedBvecsFile=*)
        	processedBvecsFile=${1#*=}
        	shift 		
        ;;  
        --processedBvalsFile=*)
        	processedBvalsFile=${1#*=}
        	shift 		
        ;;  
        --acqpFile=*)
        	acqpFile=${1#*=}
        	shift 		
        ;;  
        --indexFile=*)
        	indexFile=${1#*=}
        	shift 		
        ;;      
        --eddyVersion=*)
            eddyVersion=${1#*=}
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
        error "Unkown input argument '${1}'"
        shift 		
        ;;   			
    esac
	done	
}

#######################################
# Parse input
#######################################
computeReferenceB0(){
    
    # use b0 scan(s) to create reference volume
    for i in ${b0Scans[@]}; do
        fslroi "$dwiProcessedFile" "${dwiReferenceFile/.nii.gz/_indv_$i.nii.gz}" $((i- 1)) 1
    done

    fslmerge -t "$dwiReferenceFile" "${dwiReferenceFile/.nii.gz/_indv_}"*.nii.gz
    fslmaths "$dwiReferenceFile" -Tmean "$dwiReferenceFile"
    rm "${dwiReferenceFile/.nii.gz/_indv_}"*.nii.gz

}

#######################################
# Main function
#######################################

trap error ERR

parse_input "$@"

dwiFiles=(${dwiFiles//,/ })
fslmerge -t "$dwiProcessedFile" "${dwiFiles[@]}"

# TODO: check FSL version > 5.0.9

# Load b0-scans list and parse them into an array
b0Scans=(${b0Scans//,/ })

dwiB0ScansFile=${dwiProcessedFile/.nii.gz/_b0.nii.gz}
dwiB0UnwarpedScansFile=${dwiProcessedFile/.nii.gz/_b0_unwarped.nii.gz}
brainMaskFile=${dwiProcessedFile/.nii.gz/_mask.nii.gz} # DO NOT CHANGE: HARDCODED IN FSL

# -------------------------------------
# Topup
# -------------------------------------

# Extract b0-scan(s) and put them (temporary) in dwiReferenceFile
for i in ${b0Scans[@]}; do
	fslroi "$dwiProcessedFile" "${dwiB0ScansFile/.nii.gz/_indv_$i.nii.gz}" $((i - 1)) 1
done
fslmerge -t "$dwiB0ScansFile" "${dwiB0ScansFile/.nii.gz/_indv_}"*.nii.gz
rm "${dwiB0ScansFile/.nii.gz/_indv_}"*.nii.gz

# Run topup
topup -v --imain="$dwiB0ScansFile" --config=b02b0.cnf \
--datain="$acqpFile" --out="${dwiProcessedFile/.nii.gz/}" \
--iout="$dwiB0UnwarpedScansFile" --subsamp=1
rm "${dwiB0ScansFile/.nii.gz/}".{nii.gz,topup_log}

# -------------------------------------
# Eddy
# -------------------------------------

# Create reference b0-image from unwarped b0-scans
fslmaths "$dwiB0UnwarpedScansFile" -Tmean "$dwiReferenceFile"
bet "$dwiReferenceFile" "$dwiProcessedFile" -m -n # creates $brainMaskFile
rm "$dwiReferenceFile" "$dwiB0UnwarpedScansFile"

# Run eddy current correction
$eddyVersion -v --imain="$dwiProcessedFile" \
             --mask="$brainMaskFile" \
             --acqp="$acqpFile" \
             --index="$indexFile" \
             --bvecs="$processedBvecsFile" --bvals="$processedBvalsFile" \
             --topup="${dwiProcessedFile/.nii.gz/}" \
             --out="$dwiProcessedFile"

# clean up
# rm ${dwiProcessedFile/.nii.gz/_fieldcoef.nii.gz} ${dwiProcessedFile/.nii.gz/_movpar.txt}
rm "$brainMaskFile"

mv "${dwiProcessedFile}".eddy_rotated_bvecs "$processedBvecsFile"

# -------------------------------------
# DWI and Freesurfer registration
# -------------------------------------

# Create new eddy current and movement corrected reference b0-scan
computeReferenceB0

# Compute registration matrix from FS to subject space
bbregister --s "$freesurferDir" --mov "$dwiReferenceFile" --reg "$registrationMatrixFile" \
--dti --init-fsl
rm "`basename $registrationMatrixFile .dat`".{dat.mincost,dat.param,dat.sum,log}

# Register freesurfer segmentation to reference volume
mri_label2vol --seg "${freesurferDir}/mri/aseg.mgz" --temp "$dwiReferenceFile" \
--reg "$registrationMatrixFile" --o "$segmentationFile"

# TODO: add registration of MTR weighted image

exit 0
