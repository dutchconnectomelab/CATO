function functional_preprocessing(configParams)
% FUNCTIONAL_PREPROCESSING   Pre-process MRI data.
%
%   functional_preprocessing(CONFIGPARAMS) takes FreeSurfer and rs-fMRI data
%   and generates intermediate files for the consecutive processing steps
%   according to parameters specified in CONFIGPARAMS.
%
%   Note: CONFIGPARAMS should be the "computed" configuration parameters, with
%   hardcoded subject directory. This function is expected to be called by
%   functional_pipeline. To run this specific step use:
%   functional_pipeline(SUBJECTDIR, 'reconstructionSteps', 'functional_preprocessing')

%% Initialization

status.functional_preprocessing = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---functional_preprocessing started----\n');

%% Run preprocessing scripts

% Setup freesurfer and other stuff?
preprocessingScript = configParams.functional_preprocessing.preprocessingScript;

% Prepare input arguments
inputArguments = [struct2list(configParams.functional_preprocessing)' ...
    struct2list(configParams.general)']; % TODO: remove general

% transform logical (sliceTimingCorrection) to char.
inputArguments(2,:) = cellfun(@jsonencode, inputArguments(2,:), ...
    'UniformOutput', false); 

inputArguments = inputArguments(:, ...
    ~strcmp(inputArguments(1, :), 'preprocessingScripts'));

inputArguments = sprintf('--%s=%s \\\n\t', inputArguments{:});

% Execute preprocesing script 
cmd = ['bash ' preprocessingScript, ' ', inputArguments(1:end-3)];
exitCode = system(cmd);

if exitCode ~=0
    error('CATO:functional_preprocessing:errorInPreprocessScript', ...
        'Error in executing preprocessing script:\n %s', ...
        cmd);
end

%% Clean up

status.functional_preprocessing = 'finished';
updateStatus(configParams.general.statusFile, status);

fprintf('---functional_preprocessing finished----\n');

