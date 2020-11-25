function structural_preprocessing(configParams)
% STRUCTURAL_PREPROCESSING   Pre-process MRI data.
%
%   structural_preprocessing(CONFIGPARAMS) takes FreeSurfer and DWI data
%   and generates intermediate files for the consecutive processing steps
%   according to parameters specified in CONFIGPARAMS.
%
%   Note: CONFIGPARAMS should be the "computed" configuration parameters, with
%   hardcoded subject directory. This function is expected to be called by
%   structural_pipeline. To run this specific step use:
%   structural_pipeline(SUBJECTDIR, 'reconstructionSteps', 'structural_preprocessing')

%% Initialization

status.structural_preprocessing = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---Preprocessing started----\n');

% Initialize parameters:
BVAL_THRESHOLD = 10;

dwiFile = configParams.structural_preprocessing.dwiFile;
dwiFileReversed = configParams.structural_preprocessing.dwiReversedFile;
processedBvecsFile = configParams.structural_preprocessing.processedBvecsFile;
processedBvalsFile = configParams.structural_preprocessing.processedBvalsFile;
acqpFactor = configParams.structural_preprocessing.acqpFactor;
revPhaseEncDim = configParams.structural_preprocessing.revPhaseEncDim;
dwiB0OnlyReversed = configParams.structural_preprocessing.dwiB0ReversedFile;        
indexFile = configParams.structural_preprocessing.indexFile;
acqpFile = configParams.structural_preprocessing.acqpFile;
rawBvalsFile = configParams.structural_preprocessing.rawBvalsFile;
rawBvecsFile = configParams.structural_preprocessing.rawBvecsFile;

%% Load bvals and bvecs

% Load validate and standardize b-values
bvals = dlmread(rawBvalsFile);
validateattributes(bvals, ...
    {'numeric'}, {'nonempty', 'vector', 'nonnegative'}, ...
    mfilename, 'rawBvals');
            
bvals = bvals(:);

% Load validate and standardize b-vectors
bvecs = dlmread(rawBvecsFile);
validateattributes(bvecs, ...
    {'numeric'}, {'nonempty', '2d'}, ...
    mfilename, 'rawBvals');

assert(any(size(bvecs) == 3), ...
    'bvecs must match number of bVals N and be size Nx3 or 3xN matrix');

if size(bvecs, 1) == 3
    bvecs = bvecs';
end

assert(isequal(length(bvals), size(bvecs, 1)), ...
    'CATO:structural_pipeline:bValsAndbVecsDoNotMatch', ...
    'Number of b-values (%i) and b-vectors (%i) do not match', ...
    length(bvals), size(bvecs, 1));

%% Check dimensions of DWI file and bvals and bvecs

% Check DWI files match with bvals and bvecs.
% Check DWI files match with eachother

propdwiFile = getNiftiProperties(dwiFile);
nScans = propdwiFile.dim(4);
assert(nScans == length(bvals));

if ~isempty(dwiB0OnlyReversed)
    propdwiB0OnlyReversed = getNiftiProperties(dwiB0OnlyReversed);
    assert(isequal(propdwiFile.dim(1:3), propdwiB0OnlyReversed.dim(1:3)));
    nB0ReversedScans = propdwiB0OnlyReversed.dim(4);
end

if ~isempty(dwiFileReversed)
    propdwiFileReversed = getNiftiProperties(dwiFileReversed);
    assert(isequal(propdwiFile.dim(1:4), propdwiFileReversed.dim(1:4)));
end


%% Prepare index, acqp, bvals and bvecs for Eddy and Topup


% Create acquisition parameters matrix
acqp = zeros(2, 4);
acqp(:, 4) = acqpFactor;
acqp(1, revPhaseEncDim) = 1; acqp(2, revPhaseEncDim) = -1;

% Case 1: One set (Minimal or Eddy preprocessing)
if isempty(dwiFileReversed) && isempty(dwiB0OnlyReversed)
    
    index = ones(1, nScans);
    
end
    
% Case 2: All scans reversed (TOPUP)
if ~isempty(dwiFileReversed) && isempty(dwiB0OnlyReversed)
        
    nB0dwi = nnz(bvals < BVAL_THRESHOLD);
    index1 = cumsum(bvals < BVAL_THRESHOLD);
    
    % Scans for the first b0-scan are associated with first b0-scan.
    index1(index1 == 0) = 1;
    index2 = index1 + max(index1);

    index = [index1; index2];
    acqp = acqp([ones(nB0dwi, 1); 2*ones(nB0dwi, 1)], :); % make a row for each b0-scan
    
    bvals = [bvals; bvals];
    bvecs = [bvecs; bvecs];
    
    configParams.structural_preprocessing.dwiFiles = ...
        strcat(dwiFileReversed, ',', dwiFile);

end

% Case 3: Only b0-scans reversed (TOPUP)
if isempty(dwiFileReversed) && ~isempty(dwiB0OnlyReversed)

    % create index file
    nB0dwi = nnz(bvals < BVAL_THRESHOLD);
    nB0ReversedScans = propdwiB0OnlyReversed.dim(4);
    
    index1 = (1:nB0ReversedScans)';
    index2 = cumsum(bvals < BVAL_THRESHOLD);
    index2(index2 == 0) = 1; % Scans for the first b0-scan are associated with first b0-scan.  
    index2 = index2 + max(index1);

    index = [index1; index2];
        
    % create acquisition parameter file
    acqp = acqp([ones(nB0ReversedScans, 1); 2*ones(nB0dwi, 1)], :);
    
    % create bvals and bvecs files
    bvals1 = zeros(nB0ReversedScans, 1);
    bvecs1 = zeros(nB0ReversedScans, 3);
    
    bvals = [bvals1; bvals];
    bvecs = [bvecs1; bvecs];
    
    configParams.structural_preprocessing.dwiFiles = ...
        strcat(dwiB0OnlyReversed, ',', dwiFile);

end

if ~isempty(dwiFileReversed) && ~isempty(dwiB0OnlyReversed)
   error('CATO:structural_preprocessing:TopupUnclear', ...
       'Either dwiFileReversed or dwiB0OnlyReversed, or neither variables should be defined.');
end

% Save standardized b-values and b-vectors
dlmwrite(processedBvecsFile, bvecs, 'delimiter', ' ');
dlmwrite(processedBvalsFile, bvals, 'delimiter', ' ');
dlmwrite(acqpFile, acqp, 'delimiter', ' ');
dlmwrite(indexFile, index, 'delimiter', ' ');

if index(1) == false
   warning('CATO:structural_preprocessing', ...
       'First DWI scan is not a b0-scan, this might raise errors.'); 
end

%% Run preprocessing scripts

% Make a list of the b0-scans
b0Scans = find(bvals < BVAL_THRESHOLD)';
configParams.structural_preprocessing.b0Scans = ...
    strrep(num2str(b0Scans), '  ', ',');

% Setup freesurfer and other stuff?
preprocessingScript = configParams.structural_preprocessing.preprocessingScript;

% Prepare input arguments
% Include general for freesurferDir
inputArguments = [struct2list(configParams.structural_preprocessing)' ...
    struct2list(configParams.general)']; 
inputArguments = inputArguments(:, cellfun(@ischar, inputArguments(2,:)));
indx = ~strcmp(inputArguments(1, :), 'preprocessingScripts');
inputArguments = inputArguments(:, indx);

inputArguments = sprintf('--%s="%s" \\\n\t', inputArguments{:});

% Standard preprocessing script parameters:
%     b0Scans
%     dwiProcessedFile
%     dwiFiles
%     refB0File
%     freesurferDir
%     registrationMatrix
%     segmentationFile
%     b0Scans
%     processedBvecsFile
%     processedBvalsFile
%     acqpFile
%     indexFile
%     eddyVersion

% Execute preprocesing script 
cmd = ['"', preprocessingScript, '" ', inputArguments(1:end-3)];
exitCode = system(cmd);

if exitCode ~=0
    error('CATO:structural_preprocessing:errorInPreprocessScript', ...
        'Error in executing preprocessing script:\n %s', ...
        cmd);
end

%% Clean up

status.structural_preprocessing = 'finished';
updateStatus(configParams.general.statusFile, status);

fprintf('---Preprocessing finished----\n');

