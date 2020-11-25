function setup_environment(configParams)
% SETUP_ENVIRONMENT Setup shell environment.
%
%   setup_environment(configParams) sets FreeSurfer (FREESURFER_HOME,
%   SUBJECTS_DIR) and FSL (FSLDIR) values in the shell environment.
%   Further, it sets the Bash startup file (BASH_ENV) to the
%   setup_environment.sh script that sources the FreeSurfer
%   (SetUpFreeSurfer.sh) and FSL (fsl.sh) setup scripts.
%
%   INPUT VARIABLES
%   configParams:
%   COnfiguration parameters.

setenv('FREESURFER_HOME');
setenv('SUBJECTS_DIR');
setenv('FSLDIR');
setenv('BASH_ENV');

testFile = fullfile(configParams.general.freesurferRootDir, 'bin', 'recon-all');
if ~isfile(testFile)
    warning(['general.freesurferRootDir (''%s'') seems corrupt ', ...
        '(''%s'' is not a file).'], configParams.general.freesurferRootDir, testFile);
end

testFile = fullfile(configParams.general.fslRootDir, 'bin', 'bet');
if ~isfile(testFile)
    warning(['general.fslRootDir (''%s'') seems corrupt ', ...
        '(''%s'' is not a file).'], configParams.general.fslRootDir, testFile);
end

% source Freesurfer
if ~isempty(configParams.general.freesurferRootDir) && ...
        ~strcmpi(configParams.general.freesurferRootDir, 'none') 
    setenv('FREESURFER_HOME', configParams.general.freesurferRootDir);
    setenv('SUBJECTS_DIR', '.');
end

if ~isempty(configParams.general.fslRootDir) && ...
        ~strcmpi(configParams.general.fslRootDir, 'none') 
    setenv('FSLDIR', configParams.general.fslRootDir);
end

setenv('BASH_ENV', fullfile(fileparts(mfilename('fullpath')), ...
    'setup_environment.sh'));

end
