function [fmriPartialvol, hdr] = load_nifti_fmri(configParams, brainMask)
% LOAD_NIFTI_MRI Load voxels in mask from fmriProcessedFile.
%
%   load_nifti_fmri(configParams, brainMask) first unzips the
%   fmriProcessedFile to a temporary file fmriProcessedTMPFile. Thereafter
%   voxels with values TRUE in the brainMask are loaded into fmri.


%% If no voxels are selected in brain mask, then load only nifti header

if isempty(brainMask)
    
[~, ~, ext] = fileparts(configParams.functional_preprocessing.fmriProcessedFile);

if strcmp(ext, '.gz')
    
    % Decompress only first 350 bytes that contain the header (348 bytes).
    tmpFile = tempname(configParams.general.TEMPDIR);
    cmd = ['gunzip -c "' configParams.functional_preprocessing.fmriProcessedFile ...
        '" | head -c 350 > "' tmpFile '"'];
    [status, output] = system(cmd);
    
    if status ~= 0
        error(['Error unzipping fmriProcessedFile:\n ''' cmd '''.\n%s'], output);
    end
end

fmriPartialvol = [];
hdr = load_nifti_hdr(tmpFile);

delete(tmpFile);

return

end

%% Unzip nifti file

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
fmriPartialvol = single(fmri.partialvol);
hdr = rmfield(fmri, 'partialvol');