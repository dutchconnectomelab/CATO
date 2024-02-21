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
fprintf('---structural_preprocessing started----\n');

% Initialize parameters:
bValueScalingTol = configParams.general.bValueScalingTol;
bValueZeroThreshold = configParams.general.bValueZeroThreshold;

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
synb0File = configParams.structural_preprocessing.synb0File;

%% Load bvals and bvecs

% Load gradient table
gtab = load_gtab(rawBvalsFile, rawBvecsFile, ...
    bValueZeroThreshold, bValueScalingTol);

bvecs = gtab.bvecs;
bvals = gtab.bvals;

%% Check dimensions of DWI file and bvals and bvecs

% Check DWI files match with bvals and bvecs.
% Check DWI files match with eachother

propdwiFile = load_nifti_hdr_fast(dwiFile);
nScans = propdwiFile.dim(5);
assert(nScans == length(bvals));

if ~isempty(dwiB0OnlyReversed)
    propdwiB0OnlyReversed = load_nifti_hdr_fast(dwiB0OnlyReversed);
    assert(isequal(propdwiFile.dim(2:4), propdwiB0OnlyReversed.dim(2:4)));
    nB0ReversedScans = propdwiB0OnlyReversed.dim(5);
end

if ~isempty(dwiFileReversed)
    propdwiFileReversed = load_nifti_hdr_fast(dwiFileReversed);
    assert(isequal(propdwiFile.dim(2:5), propdwiFileReversed.dim(2:5)));
end

% When using synboFile
% could be expanded to every topup run
if ~isempty(synb0File)
    propsynb0File = load_nifti_hdr_fast(synb0File);
    % check if dimensons have an even number to enable topup 
    if any(mod(propsynb0File.dim(2:4), 2)) || any(mod(propdwiFile.dim(2:4), 2))
        error('CATO:structural_preprocessing:Synb0Dimensions', ...
            'synb0File and dwiFile should have even dimensions to run topup.');
    end
end


%% Prepare index, acqp, bvals and bvecs for Eddy and Topup
% - deafult behavior is to generate acqp and index files
% - enable user defined acqp and index if synb0File is defined

% Set acqp based on revPhaseEncDim and acqpFactor (default)
acqp = zeros(2, 4);
acqp(:, 4) = acqpFactor;
acqp(1, revPhaseEncDim) = 1; acqp(2, revPhaseEncDim) = -1;


% Case 1: One set (Minimal or Eddy preprocessing)
if isempty(dwiFileReversed) && isempty(dwiB0OnlyReversed) && isempty(synb0File)
    
    index = ones(1, nScans);
    
end
    
% Case 2: All scans reversed (TOPUP)
if ~isempty(dwiFileReversed) && isempty(dwiB0OnlyReversed)
        
    nB0dwi = nnz(bvals == 0);
    index1 = cumsum(bvals == 0);
    
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
    nB0dwi = nnz(bvals == 0);
    
    index1 = (1:nB0ReversedScans)';
    index2 = cumsum(bvals == 0);
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

% Case 4: SynB0-DISCO output is available to run TOPUP and eddy
if ~isempty(synb0File) && isempty(dwiFileReversed) && isempty(dwiB0OnlyReversed)
    % check if index and acqp files are defined
    if exist(acqpFile, 'file') && exist(indexFile, 'file')
        acqp = dlmread(acqpFile);
        index = dlmread(indexFile);

        % check if acqp has at least 1 zero in the 4th column
        if ~any(acqp(:, 4) == 0)
            error('CATO:structural_preprocessing:TopupUnclear', ...
                'synb0File was set. acqpFile should have at least one echotime of zero (4th column).');
        end

        % if outputDir does not contain acqpFile and indexFile
        % copy them to outputDir to keep track of the used files
        if ~strcmp(configParams.general.outputDir, fileparts(acqpFile))
            copyfile(acqpFile, fullfile(configParams.general.outputDir, 'acqp.txt'));
            copyfile(indexFile, fullfile(configParams.general.outputDir, 'index.txt'));
        end

    elseif ~exist(acqpFile, 'file') || ~exist(indexFile, 'file')
        error('CATO:structural_preprocessing:TopupUnclear', ...
            'acqpFile and indexFile should be defined if synb0File is defined.');
    end



% error if conflicting inputs (synb0File defined, but also dwiFileReversed or dwiB0OnlyReversed)
elseif ~isempty(synb0File) && (~isempty(dwiFileReversed) || ~isempty(dwiB0OnlyReversed))
    error('CATO:structural_preprocessing:TopupUnclear', ...
       'Either synb0File or dwiFileReversed or dwiB0OnlyReversed, or neither variables should be defined.');
end


% error if conflicting inputs (dwiFileReversed and dwiB0OnlyReversed)
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
b0Scans = find(bvals == 0)';
configParams.structural_preprocessing.b0Scans = ...
    strrep(num2str(b0Scans), '  ', ',');

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
%     synb0File
%     b0Scans
%     processedBvecsFile
%     processedBvalsFile
%     acqpFile
%     indexFile
%     eddyVersion

% Execute preprocesing script
fprintf('Execute preprocessing script:\n%s\n', preprocessingScript);
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

fprintf('---structural_preprocessing finished----\n');

