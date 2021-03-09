function reconstruction_diffusion(configParams)
% RECONSTRUCTION_DIFFUSION   Estimate diffusion peaks and diffusion measures.
%
%   reconstruction_diffusion(CONFIGPARAMS) reconstruct diffusion peaks
%   based on diffusion MRI data based on the parameters specified in
%   CONFIGPARAMS.

%% Initialization

status.reconstruction_diffusion = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_diffusion started----\n');

reconMethods = configParams.general.reconstructionMethods;

dwiProcessedFile = configParams.structural_preprocessing.dwiProcessedFile;
dwiReferenceFile = configParams.structural_preprocessing.dwiReferenceFile;
segmentationFile = configParams.structural_preprocessing.segmentationFile;
processedBvalsFile = configParams.structural_preprocessing.processedBvalsFile;
processedBvecsFile = configParams.structural_preprocessing.processedBvecsFile;

bValueScalingTol = configParams.general.bValueScalingTol;
bValueZeroThreshold = configParams.general.bValueZeroThreshold;

diffusionPeaksFile = configParams.reconstruction_diffusion.diffusionPeaksFile;
diffusionMeasuresFile = configParams.reconstruction_diffusion.diffusionMeasuresFile;

nonlinearitiesFlag = configParams.reconstruction_diffusion.gradientNonlinearities.correctNonlinearities;
nonlinearitiesFile = configParams.reconstruction_diffusion.gradientNonlinearities.nonlinearitiesFile;

% Load processed DWI file.
try
    dwi = load_nifti(dwiProcessedFile);
    signalIntensities = single(dwi.vol);
    dim = size(dwi.vol);
    clear dwi
catch
    error('Cannot load dwiProcessedFile (%s).', dwiProcessedFile);
end
signalIntensities = reshape(signalIntensities, [], size(signalIntensities, 4));

% Load segmentation file.
try
    segmentation = load_nifti(segmentationFile);
    maskBrain = segmentation.vol > 0;
catch
    error('Cannot load segmentationFile (%s).', dwiProcessedFile);
end

assert(isequal(size(maskBrain), dim(1:3)), ...
    'segmentationFile (%s) and dwiProcessedFile (%s) differ in voxel dimensions', ...
    segmentationFile, dwiProcessedFile);

% Load gradient table
gtab = load_gtab(processedBvalsFile, processedBvecsFile, bValueZeroThreshold, bValueScalingTol);

% Import voxel wise gradient correction.
if nonlinearitiesFlag
    try
        nonlinearities = load_nifti(nonlinearitiesFile);
        nonlinearities = nonlinearities.vol;
    catch
        error('Cannot load nonlinearitiesFile (%s).', nonlinearitiesFile);
    end
    
    nonlinearities = reshape(nonlinearities, [], size(nonlinearities, 4));
end

%% Run diffusion reconstruction algorithms.

% Select voxels
indxV = find(maskBrain(:) & (mean(signalIntensities > 0, 2) >= 0.9));
nVoxels = size(signalIntensities, 1);

%% DTI reconstruction
if any(strcmpi(reconMethods, 'dti'))
    disp('DTI reconstruction started');
    
    [thresCondNum, thresVarProjScores] = thresholdAssistant(gtab);
    
    diffusionTensors = nan(nVoxels, 6);
    if ~nonlinearitiesFlag
        diffusionTensors(indxV, :) = iRESTORE(signalIntensities(indxV, :), gtab, ...
            thresCondNum, thresVarProjScores);
    else
        diffusionTensors(indxV, :) = iRESTORE(signalIntensities(indxV, :), gtab, ...
            thresCondNum, thresVarProjScores, nonlinearities);
    end
    
    diffusionPeaks = tensor2peak(diffusionTensors); %#ok
    [diffusionMeasures, weightDescriptions] = tensor2measure(diffusionTensors); %#ok
    diffusionMeasures = reshape(diffusionMeasures, [dim(1:3) 4]); %#ok
    
    save(strrep(diffusionPeaksFile, 'METHOD', 'dti'), 'diffusionPeaks');
    save(diffusionMeasuresFile, 'diffusionMeasures', 'weightDescriptions');
    
    clear diffusionPeaks diffusionMeasures
    
    disp('DTI reconstruction finished');
end

%% CSD reconstruction
if any(strcmpi(reconMethods, 'csd'))
    
    disp('CSD reconstruction started');
        
    CCRegions = configParams.reconstruction_diffusion.CSD.CCRegions;
    FAThreshold = configParams.reconstruction_diffusion.CSD.FAThreshold;
    
    lambda = configParams.reconstruction_diffusion.CSD.lambda;
    tau = configParams.reconstruction_diffusion.CSD.tau;
    shOrder = configParams.reconstruction_diffusion.CSD.shOrder;
    
    outputPeaks = configParams.reconstruction_diffusion.CSD.outputPeaks;
    minPeakRatio = configParams.reconstruction_diffusion.CSD.minPeakRatio;
    maxPeaks = configParams.reconstruction_diffusion.CSD.maxPeaks;
    
    CCMask = ismember(segmentation.vol, CCRegions);
    if ~any(CCMask)
        error('CATO:reconstruction_diffusion:noCCMaskVoxels', ...
            ['No voxel with segmentation codes matching the Corpus Callosum ', ...
            'CCRegions: %s'], mat2str(CCRegions));
    end
    
    if isfile(diffusionMeasuresFile)
        
        try
            data = load(diffusionMeasuresFile);
            diffusionMeasures = data.diffusionMeasures;
            weightDescriptions = data.weightDescriptions;
            clear data
        catch
            error(['Cannot load diffusionMeasuresFile (including variables', ...
                ' diffusionMeasures and weightDescriptions) (%s).'], diffusionMeasuresFile);
        end
        
        assert(all(ismember({'fractional anisotropy', ...
            'radial diffusivity', ...
            'axial diffusivity'}, weightDescriptions)), ...
            'CATO:DTIMeasuresNotFound', ...
            'DTI measures not found in diffusionMeasuresFile. Reconstruct DWI signal using DTI prior to running CSD');
    else
        error('CATO:diffusionMeasuresFileNotFound', ...
            'Diffusion measures file not found. Reconstruct DWI signal using DTI prior to running CSD');
    end
    
    response = auto_response(diffusionMeasures, weightDescriptions, signalIntensities, ...
        gtab.bvals, CCMask, FAThreshold);
    
    diffusionPeaks = zeros(nVoxels, 3, outputPeaks);
    
    if ~nonlinearitiesFlag
        diffusionPeaks(indxV, :, :) = srcsd(signalIntensities(indxV, :), gtab, response, ...
            shOrder, lambda, tau, outputPeaks, minPeakRatio, maxPeaks); %#ok
    else
        diffusionPeaks(indxV, :, :) = srcsd(signalIntensities(indxV, :), gtab, response, ...
            shOrder, lambda, tau, outputPeaks, minPeakRatio, maxPeaks, nonlinearities); %#ok
    end
    
    save(strrep(diffusionPeaksFile, 'METHOD', 'csd'), 'diffusionPeaks');
    
    clear diffusionPeaks
    
    disp('CSD reconstruction finished');
    
end

%% GQI reconstruction
if any(strcmpi(reconMethods, 'gqi'))
    
    disp('GQI reconstruction started');
    
    outputPeaks = configParams.reconstruction_diffusion.GQI.outputPeaks;
    minPeakRatio = configParams.reconstruction_diffusion.GQI.minPeakRatio;
    maxPeaks = configParams.reconstruction_diffusion.GQI.maxPeaks;
    GQIratio = configParams.reconstruction_diffusion.GQI.meanDiffusionDistanceRatio;
    
    diffusionPeaks = zeros(nVoxels, 3, outputPeaks);
    GFA = nan(nVoxels, 1);
    QA = nan(nVoxels, outputPeaks);
    
    if ~nonlinearitiesFlag
        [diffusionPeaks(indxV, :, :) , GFA(indxV) , QA(indxV, :)] = gqi(signalIntensities(indxV, :), gtab, ...
            GQIratio, outputPeaks, minPeakRatio, maxPeaks); %#ok
    else
        [diffusionPeaks(indxV, :, :) , GFA(indxV) , QA(indxV, :)] = gqi(signalIntensities(indxV, :), gtab, ...
            GQIratio, outputPeaks, minPeakRatio, maxPeaks, nonlinearities); %#ok
    end
    
    save(strrep(diffusionPeaksFile, 'METHOD', 'gqi'), 'diffusionPeaks');
    
    NGFA = GFA ./ nanmax(GFA); % normalized GFA
    
    % Save GQI diffusion measures
    diffusionMeasuresGQI = cat(2, GFA, NGFA, QA);
    diffusionMeasuresGQI = reshape(diffusionMeasuresGQI, ...
        [dim(1:3) size(diffusionMeasuresGQI, 2)]);
    weightDescriptionsGQI = [{'generalized fractional anisotropy'}; ...
        {'normalized generalized fractional anisotropy'}; ...
        compose('quantitative anisotropy peak #%i', 1:outputPeaks)'];
    
    if isfile(diffusionMeasuresFile)
        
        try
            data = load(diffusionMeasuresFile);
            diffusionMeasures = data.diffusionMeasures;
            weightDescriptions = data.weightDescriptions;
        catch
            error(['Cannot load diffusionMeasuresFile (including variables', ...
                'diffusionMeasures and weightDescriptions) (%s).'], diffusionMeasuresFile);
        end
        
        matchDimensions = size(diffusionMeasures) == size(diffusionMeasuresGQI);
        assert(all(matchDimensions(1:3)), ...
            ['Dimensions diffusionMeasuresFile (%s) do not match ', ...
            'reconstructed GQI diffusion measures dimensions'], diffusionMeasuresFile);
        
        % delete old GFA and QA measures (if exist)
        indx = ismember(weightDescriptions, weightDescriptionsGQI);
        diffusionMeasures = diffusionMeasures(:, :, :, ~indx);
        weightDescriptions = weightDescriptions(~indx);
        
    else
        
        diffusionMeasures = [];
        weightDescriptions = {};
        
    end
    
    diffusionMeasures = cat(4, diffusionMeasures, diffusionMeasuresGQI); %#ok
    weightDescriptions = cat(1, weightDescriptions, weightDescriptionsGQI); %#ok
    
    save(diffusionMeasuresFile, 'diffusionMeasures', 'weightDescriptions');
    
    clear diffusionMeasures diffusionMeasuresGQI diffusionPeaks
    
    disp('GQI reconstruction finished');
    
end

%% GQI_DTI reconstruction
if any(strcmpi(reconMethods, 'gqi_dti'))
    
    disp('GQI_DTI reconstruction started');
    
    try
        diffusionPeaksDTI = load(strrep(diffusionPeaksFile, 'METHOD', 'dti'), 'diffusionPeaks');
        diffusionPeaksDTI = diffusionPeaksDTI.diffusionPeaks;
    catch
        error('CATO:CannotOpenDTIPeaks', ...
            ['DTI peaks file cannot be opened. Reconstruct DWI signal ', ...
            'using DTI prior to running GQI_DTI']);
    end
    
    try
        diffusionPeaksGQI = load(strrep(diffusionPeaksFile, 'METHOD', 'gqi'), 'diffusionPeaks');
        diffusionPeaksGQI = diffusionPeaksGQI.diffusionPeaks;
    catch
        error('CATO:CannotOpenGQIPeaks', ...
            ['GQI peaks file cannot be opened. Reconstruct DWI signal ', ...
            'using GQI prior to running GQI_DTI']);
    end
    
    indxCrossingFibers = sum(any(diffusionPeaksGQI, 2),3) > 1;
    
    diffusionPeaks = zeros(size(diffusionPeaksGQI));
    diffusionPeaks(~indxCrossingFibers, :, 1) = diffusionPeaksDTI(~indxCrossingFibers, :);
    diffusionPeaks(indxCrossingFibers, :, :) = diffusionPeaksGQI(indxCrossingFibers, :, :); %#ok
    
    save(strrep(diffusionPeaksFile, 'METHOD', 'gqi_dti'), 'diffusionPeaks');
    
    disp('GQI_DTI reconstruction finished');
    
end

%% CSD_DTI reconstruction
if any(strcmpi(reconMethods, 'csd_dti'))
    
    disp('CSD_DTI reconstruction started');
    try
        diffusionPeaksDTI = load(strrep(diffusionPeaksFile, 'METHOD', 'dti'), 'diffusionPeaks');
        diffusionPeaksDTI = diffusionPeaksDTI.diffusionPeaks;
    catch
        error('CATO:CannotOpenDtiPeaks', ...
            'DTI peaks file cannot be opened. Reconstruct DWI signal using DTI prior to running CSD_DTI');
    end
    
    try
        diffusionPeaksCSD = load(strrep(diffusionPeaksFile, 'METHOD', 'csd'), 'diffusionPeaks');
        diffusionPeaksCSD = diffusionPeaksCSD.diffusionPeaks;
    catch
        error('CATO:CannotOpenCSDPeaks', ...
            'CSD peaks file cannot be opened. Reconstruct DWI signal using CSD prior to running CSD_DTI');
    end
    
    indxCrossingFibers = sum(any(diffusionPeaksCSD, 2),3) > 1;
    
    diffusionPeaks = zeros(size(diffusionPeaksCSD));
    diffusionPeaks(~indxCrossingFibers, :, 1) = diffusionPeaksDTI(~indxCrossingFibers, :);
    diffusionPeaks(indxCrossingFibers, :, :) = diffusionPeaksCSD(indxCrossingFibers, :, :); %#ok
    
    save(strrep(diffusionPeaksFile, 'METHOD', 'csd_dti'), 'diffusionPeaks');
    
    clear diffusionPeaksDTI diffusionPeaksCSD diffusionPeaks
    
    disp('CSD_DTI reconstruction finished');
    
end

%% Export to NIFTI

if configParams.reconstruction_diffusion.exportNifti.exportNifti
    
    exportMeasures = configParams.reconstruction_diffusion.exportNifti.measures;
    exportNiftiFile = configParams.reconstruction_diffusion.exportNifti.diffusionMeasuresFileNifti;
    
    try
        data = load(diffusionMeasuresFile);
        diffusionMeasures = data.diffusionMeasures;
        weightDescriptions = data.weightDescriptions;
    catch
        error(['Cannot load diffusionMeasuresFile (including variables', ...
            ' diffusionMeasures and weightDescriptions) (%s).'], diffusionMeasuresFile);
    end
    
    diffusionMeasuresNifti = load_nifti(dwiReferenceFile);
    
    for iM = 1:length(exportMeasures)
        [~, indxMeasures] = ismember(exportMeasures{iM}, weightDescriptions);
        if indxMeasures == 0
            error(['Diffusion measure %s not available in ', ...
                'diffusionmeasuresFile (%s). Make sure you run all ', ...
                'reconstruction methods (e.g. DTI) necessary.'], ...
                exportMeasures{iM}, diffusionMeasuresFile);
        end
        diffusionMeasuresNifti.vol = diffusionMeasures(:, :, :, indxMeasures);
        save_nifti(diffusionMeasuresNifti, strrep(exportNiftiFile, ...
            'MEASURE', genvarname(strrep(exportMeasures{iM}, ' ', '_'))));
    end
    
    fprintf('Diffusion measures (%s) exported to Nifti files.\n', ...
        strjoin(exportMeasures, ', '));
    
    
end

%% Finish

status.reconstruction_diffusion = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_diffusion finished----\n');


