#!/bin/bash

# Setup FSL
if [ -n "$FSLDIR" ];
then
	PATH=${FSLDIR}/bin:${PATH}
	export PATH
	export FSL_DIR=$FSLDIR
	source ${FSLDIR}/etc/fslconf/fsl.sh
fi

# Setup Freesurfer
if [ -n "$FREESURFER_HOME" ];
then
	source $FREESURFER_HOME/SetUpFreeSurfer.sh
fi
