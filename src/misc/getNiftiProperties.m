function props = getNiftiProperties(niftiFile)
% GETNIFTIPROPERTIES Read properties of a Nifti-file.
%
%   props = getNiftiProperties(niftiFile) reads, quickly, the dimensions,
%   size of the dimensions, and orientation from niftiFile using mri_info.

% check mri_info is added to path.
[status, output] = system('command -v mri_info');
if status ~= 0
    error(['Cannot find mri_info. Make sure FSL and Freesurfer ', ...
        'directories are correct and use setup_environment']);
end

% Dimensions
cmd = ['mri_info --dim "' niftiFile '" 2>/dev/null | tail -n1'];
[status, output] = system(cmd);
if status ~= 0
        error(['Error in:\n''' cmd '''.\n%s'], output);

end

try
    props.dim = str2num(output);
    assert(isvector(props.dim));
catch
    error('Nifti header cannot be parsed. Unexpected output mri_info:\n%s', output);
end

% Size dimensions
cmd = ['mri_info --res "' niftiFile '" 2>/dev/null | tail -n1'];
[status, output] = system(cmd);
if status ~= 0
        error(['Error in:\n''' cmd '''.\n%s'], output);

end
try
    props.res = str2num(output);
    assert(isvector(props.res));
catch
    error('Nifti header cannot be parsed. Unexpected output mri_info:\n%s', output);
end

% Orientation
cmd = ['mri_info --orientation "' niftiFile '" 2>/dev/null | tail -n1'];
[status, output] = system(cmd);
if status ~= 0
        error(['Error in:\n''' cmd '''.\n%s'], output);
end

try
    props.orientation = strip(output);
    assert(length(props.orientation) == 3);
catch
    error('Nifti header cannot be parsed. Orientation (%s) incorrect', output);
end

end
