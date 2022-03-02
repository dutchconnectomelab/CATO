function configParams = structural_pipeline(subjectDir, varargin)
% STRUCTURAL_PIPELINE   Reconstruct structural connectivity.
%
%   structural_pipeline(SUBJECTDIR) reconstructs structural connectivity of
%   a subjects T1 and DWI data in the SUBJECTDIR directory using default
%   parameters.
%
%   structural_pipeline(SUBJECTDIR, 'Param1', VAL1, 'Param2', VAL2, ...)
%   specifies additional parameter name-value pairs chosen from:
%
%       'reconstructionSteps'   The reconstruction steps that should be
%                               executed. RECONSTRUCTIONSTEPS is a single 
%                               step or a cell array containing one or more
%                               steps.
%                               Default: all pipeline steps
%       'configurationFile'     File with reconstruction parameters.
%                               Default: 'none'
%       'runType'               Behavior when re-running CATO. Possible 
%                               values are:
%                               'none'      - Execute only steps that
%                                             have not been started earlier.
%                               'continue'  - Execute only steps that
%                                             have not been started earlier
%                                             or have ran or encountered
%                                             errors.
%                               'overwrite' - Execute all steps including
%                                             steps that finished successfully.
%                               Default: 'none'
%
%       Other parameters name-value pairs are as defined in the 
%       configuration file. For example:
%           'general.freesurferDir', 'FreeSurfer'
%           'general.templates', {'aparc', 'lausanne120', 'lausanne250'}
%           'structural_preprocessing.rawBvalsFile', 'DWI/dwi.bvals'
%
%   For more information, see <a href="matlab:
%   web('http://www.dutchconnectomelab.nl')">the CATO documentation website</a>.
            
%% Initialization

fid = fopen('VERSION', 'r');
if fid == -1 
    catoVersion = 'unknown';
else
    catoVersion = fgetl(fid);
    fclose(fid);
end

if verLessThan('matlab','9.3')
    warning('CATO:structuralPipeline:NotCompatible', ...
        ['CATO was developed and tested using MATLAB 2017b. Using an ' ...
        'older version might result in unexpected errors.']);
end

% CATO uses statistics and signal toolbox.
license('checkout','statistics_toolbox');
license('checkout','signal_toolbox');

% PARSE INPUT
[subjectDir, varargin{:}] = convertStringsToChars(subjectDir, varargin{:});

assert(ischar(subjectDir), ...
    'CATO:structural_pipeline:subjectDirNotText', ...
    'subjectDir must be a row vector of characters or string scalar.');
assert(isdir(subjectDir), ...
    'CATO:structural_pipeline:subjectDirNotDir', ...
    'subjectDir (%s) is not a directory', subjectDir);

[configFile, runType, configParamsCl] = parseVarargin(varargin{:});

%% Setup

% Setup path (note this is done before log files or other stuff).
oldPath = pwd;

try
    cd(subjectDir);
catch
    error('CATO:structural_pipeline:subjectDirNotAccessible', ...
        'Cannot open subjectDir (%s) (Name is nonexistent or not a directory).', subjectDir);
end

cleanupPath = onCleanup(@() cd(oldPath));

% To ensure different temporary file names.
rng('shuffle');

%% Combine configuration parameters

% 1. Load default configuration file
assert(exist('config_structural_default.json', 'file') == 2, ...
    ['Default configuration file (config_structural_default.json) is ', ...
    'not the search path. Make sure CATO is installed correctly.']);
configParams = readConfigFile('config_structural_default.json');

% 2. Update with user configuration file
if ~isempty(configFile)
    configFileParams = readConfigFile(configFile);
    configParams = updateParameters(configParams, configFileParams);
    configParams.general.CONFIGDIR = fileparts(configFile);
end

%  3. Update with command line parameters
configParams = updateParameters(configParams, configParamsCl);

% Fill-in variables in configuration file
[~, subjectName] = fileparts(strip(subjectDir, 'right', filesep));
assert(ischar(subjectName) & ~isempty(subjectName), ...
    'Cannot extract subject name from subject directory (%s).', subjectDir);
configParams.general.SUBJECT = subjectName;

configParams = parseConfigParams(configParams);

% Get the pipeline steps
reconStepNames = configParams.general.reconstructionSteps;

%% Initialize status file

% Select the applicable status file (.running, .error or .finished).
statusFileOld = {strrep(configParams.general.statusFile, 'STATUS', 'running'), ...
    strrep(configParams.general.statusFile, 'STATUS', 'error'), ...
    strrep(configParams.general.statusFile, 'STATUS', 'finished')};
statusFileOld = statusFileOld(isfile(statusFileOld));

if length(statusFileOld) == 1
    switch runType
        case 'overwrite'
            try
                delete(statusFileOld{1});
            catch
                error('CATO:structural_pipeline:cannotDeleteStatusFileOld', ...
                    'Cannot delete old status file (%s).', statusFileOld{1});
            end
            try
                delete(configParams.general.logFile);
            catch
                error('CATO:structural_pipeline:cannotDeleteLogFile', ...
                    'Cannot delete old log file (%s).', ...
                    configParams.general.logFile);
            end
            status = struct();
        otherwise
            status = readConfigFile(statusFileOld{1});
    end
elseif length(statusFileOld) > 1
    error('CATO:structural_pipeline:MultipleStatusfiles', ...
        ['Multiple status files exist (%s).', ...
        'Remove duplicates to ensure that pipeline runs correctly.'], ...
        strjoin(statusFileOld));
else
    status = struct();
end

%% Determine which pipeline steps should execute

% Figure out which steps should be executed.
runSteps = false(length(reconStepNames), 1);
for i = 1:length(reconStepNames)
    thisStep = reconStepNames{i};
    
    % if the status-file does not mention thisStep. Then run for all
    % runTypes.
    if ~isfield(status, thisStep)
        runSteps(i) = true;
        
        % if the status-file mentions thisStep, then figure out what to do:
    else
        switch status.(thisStep)
            case {'running', 'error'}
                if strcmp(runType, 'overwrite') || strcmp(runType, 'continue')
                    runSteps(i) = true;
                else
                    error('CATO:structural_pipeline:subjectAlreadyProcessed', ...
                        ['One or more pipeline steps have already been ', ...
                        'unsucesfully processed. Use runType ''overwrite'' ', ...
                        'or ''continue''.']);
                end
            case 'finished'
                if strcmp(runType, 'overwrite')
                    runSteps(i) = true;
                elseif strcmp(runType, 'continue')
                    runSteps(i) = false;
                else
                    error('CATO:structural_pipeline:subjectAlreadyProcessed', ...
                        ['One or more pipelinesteps have already been ', ...
                        'completed. Use runType overwrite']);
                end
            otherwise
                error('CATO:structural_pipeline:unkownStatus', ...
                    'Unknown status (%s) in statusFile (%s)', ...
                    status.(thisStep), statusFileOld{1});
        end
    end
end

reconStepNames = reconStepNames(runSteps);
reconStepNames = reconStepNames(:);

%% Setup pipeline (outputDir, logfile, #threads and welcome message)

[succes, message] = mkdir(configParams.general.outputDir);
if ~succes
    error('CATO:structural_pipeline:cannotMakeOutputDir', ...
        'Cannot make output directory outputDir (%s).\n%s', ...
        configParams.general.outputDir, message);
end

try diary(configParams.general.logFile)
catch
    error('CATO:structuralPipeline:CouldNotOpenLogFile', ...
        'Cannot not write to log file (%s)', configParams.general.logFile)
end
cleanupLog = onCleanup(@() diary('off'));

% Set maximum number of cores
if configParams.general.maxNumberCompThreads == 0
maxNumCompThreads('automatic');
else
    maxNumCompThreads(configParams.general.maxNumberCompThreads);
end
nNumCompThreads = maxNumCompThreads;
    
status.general = 'running';
updateStatus(configParams.general.statusFile, status);

timeStart = datetime('now');

fprintf('---CATO----\n');
fprintf('Version:             %s\n', catoVersion);
fprintf('Started at:          %s\n', timeStart);
fprintf('Number of threads:   %i\n', nNumCompThreads);
fprintf('Subject name:        %s\n', subjectName);
fprintf('Subject directory:   %s\n\n', subjectDir);

fprintf('Run pipeline steps:\n');
fprintf('%s\n', reconStepNames{:});
fprintf('\n');

%% Check configuration parameters (only for executed steps)

% Find errors and warnings
configErrWarnList = checkParams(configParams, [reconStepNames; {'general'}], ...
    configParams.general.parameterPropertiesFile);

% Display errors
indx = find(~cellfun(@isempty, configErrWarnList(:, 2)));
if ~isempty(indx)
    errorList = configErrWarnList(indx, [1 4 2])';
    errorList = sprintf('Parameter %s (''%s''): %s\n', errorList{:});
    error('CATO:structural_pipeline:configParamsError', ...
        'One or more parameters are incorrect:\n%s', errorList);
end

% Display warnings
indx = find(~cellfun(@isempty, configErrWarnList(:, 3)));
if ~isempty(indx)
    warningList = configErrWarnList(indx, [1 4 3])';
    warningList = sprintf('Parameter %s (''%s''): %s\n', warningList{:});
    warning('CATO:structural_pipeline:configParamsWarning', ...
        'One or more parameters are incorrect:\n%s', warningList);
end

%% Update computedConfigFile
% if a computedConfigFile already exist then update steps that are executed
% In the computed config file we save the real (i.e. parsed) variable names.

if isfile(configParams.general.computedConfigFile)
    computedConfigParams = readConfigFile(configParams.general.computedConfigFile);
    computedConfigParams = parseConfigParams(computedConfigParams);
else
    computedConfigParams = struct();
end

computedConfigParams = updateParameters(computedConfigParams, configParams);
saveConfigFile(configParams.general.computedConfigFile, computedConfigParams);

% Load the configuration file again to make sure that all substitute
% variables are replaced.
configParams = readConfigFile(configParams.general.computedConfigFile);
configParams = parseConfigParams(configParams);

%% Run pipeline steps

setup_environment(configParams)
status.general = 'finished';

for i = 1:size(reconStepNames, 1)
    try
        feval(reconStepNames{i}, configParams);
    catch ME
        MEStep = MException('CATO:structural_pipeline', ...
            'Error in %s step.', reconStepNames{i});
        ME = addCause(MEStep, ME);
        status.(reconStepNames{i}) = 'error';
        status.general = 'error';
        break
    end
end

%% Finish gracefully
updateStatus(configParams.general.statusFile, status);

switch status.general
    case 'error'
        MEcause = textwrap({ME.cause{1}.message}, 50);
        MEfile = [{ME.cause{1}.stack.name}; {ME.cause{1}.stack.line}];
        fprintf('%s\n\n', ME.message)
        fprintf('Details:\n');
        fprintf('\t%s\n', MEcause{:});
        fprintf('\tError in: %s (Line %i)\n', MEfile{:});
    case 'finished'
        timeEnd = datetime('now');
        fprintf('\n-----------------------------------\n');
        fprintf('Started at:\t%s\n', timeStart);
        fprintf('Ended at:\t%s\n', timeEnd);
        fprintf('Duration:\t%s\n', timeEnd - timeStart);
        fprintf('CATO finished without errors. \n');
end

