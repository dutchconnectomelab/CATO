function hdr = load_nifti_hdr_fast(niftiFile)
    
[~, ~, ext] = fileparts(niftiFile);

if strcmp(ext, '.gz')
    
    % Decompress only first 350 bytes that contain the header (348 bytes).
    tmpFile = tempname;
    cmd = ['gunzip -c "' niftiFile ...
        '" | head -c 350 > "' tmpFile '"'];
    [status, output] = system(cmd);
    
    if status ~= 0
        error(['Error unzipping fmriProcessedFile:\n ''' cmd '''.\n%s'], output);
    end
end

hdr = load_nifti_hdr(tmpFile);

delete(tmpFile);
