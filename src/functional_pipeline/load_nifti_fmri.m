function fmri = load_nifti_fmri(configParams, brainMask)
% LOAD_NIFTI_MRI Load voxels in mask from fmriProcessedFile.
%
%   load_nifti_fmri(configParams, brainMask) first unzips the
%   fmriProcessedFile to a temporary file fmriProcessedTMPFile. Thereafter
%   voxels with values TRUE in the brainMask are loaded into fmri.

%% Unzip
% Unzip fmriProcessedFile only if fmriProcessedTMPFile does not exist yet.
[~, ~, ext] = fileparts(configParams.functional_preprocessing.fmriProcessedFile);
if strcmp(ext, '.gz')
    fmriProcessedTMPFile = configParams.general.fmriProcessedTMPFile;
    
    if ~isfile(fmriProcessedTMPFile)
        cmd = ['gunzip -c "' configParams.functional_preprocessing.fmriProcessedFile ...
            '" > "' fmriProcessedTMPFile '"'];
        [status, output] = system(cmd);
        if status ~= 0
            error(['Error unzipping fmriProcessedFile:\n ''' cmd '''.\n%s'], output);
        end
    end
    
else
    fmriProcessedTMPFile = configParams.functional_preprocessing.fmriProcessedFile;
end

%% Read nifti partially
% Load only voxels that are selected in brainMask.
fmri = load_nifti_partially(fmriProcessedTMPFile, brainMask);
fmri = fmri.partialvol;