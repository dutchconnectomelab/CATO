function reconstruction_functional_network(configParams)
% RECONSTRUCTION_FUNCTIONAL_NETWORK Reconstruct functional connectivity.
%
%   reconstruction_functional_network(CONFIGPARAMS) reconstruct the
%   functional connectivity matrices for each template and each method
%   according to parameters specified in CONFIGPARAMS.

%% Initialization
status.reconstruction_functional_network = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_functional_network started----\n');

% Run preprocessing for several settings.
methods = fieldnames(configParams);
methods = methods(contains(methods, 'reconstruction_functional_network'));

pipelineV25Flag = false; % option for compatibility with old pipeline version.

for iMethod = 1:length(methods)
    
    %% Initializate parameters
    
    % General parameters
    thisMethodDescription = configParams.(methods{iMethod}).methodDescription;
    fmriProcessedFile = configParams.functional_preprocessing.fmriProcessedFile;
    
    segmentationFile = configParams.functional_preprocessing.segmentationFile;
    ROIsFile = configParams.general.ROIsFile;
    
    parcellationFile = configParams.parcellation.parcellationFile;
    templates = configParams.general.templates;
    templatesDir = configParams.general.templatesDir;
    
    % Regression parameters
    regression = configParams.(methods{iMethod}).regression;
    GSRFlag = regression.globalMeanRegression; % global signal regression flag
    regressionMask = regression.regressionMask;
    nMaskRegressors = length(regressionMask);
    motionParametersFile = configParams.functional_preprocessing.motionParametersFile;
    
    % Bandpass filter parameters
    minRepetitionTime = configParams.(methods{iMethod}).minRepetitionTime;
    bandpass_filter = configParams.(methods{iMethod}).bandpass_filter;
    
    % Scrubbing parameters
    scrubbing = configParams.(methods{iMethod}).scrubbing;
    motionMetricsFile = configParams.compute_motion_metrics.motionMetricsFile;
    
    % Time series parameters
    saveTimeSeriesFlag = configParams.(methods{iMethod}).saveTimeSeries;
    timeSeriesFile = configParams.(methods{iMethod}).timeSeriesFile;
    
    % Connectivity matrix parameters
    connectivityMatrixFile = configParams.(methods{iMethod}).connectivityMatrixFile;

    fprintf('method description: %s\n', thisMethodDescription);
    
    %% Prepare rs-fMRI data
    
    % Select voxels for processing that aree within the brain and have
    % signal in more than 90% of the timepoints
    segmentation = load_nifti(segmentationFile);
    segmentation = segmentation.vol;
    brainMask = segmentation(:) > 0;
    
    signalIntensities = load_nifti_fmri(configParams, brainMask);
    nTimePoints = size(signalIntensities, 2);
    
    prevalenceMask = mean(signalIntensities ~= 0, 2) >= 0.9;
    
    if ~pipelineV25Flag
    signalIntensities = signalIntensities(prevalenceMask, :);
    brainPrevalenceMask = false(size(brainMask));
    brainPrevalenceMask(brainMask) = prevalenceMask;
    else
        brainPrevalenceMask = brainMask;
    end
    
    %% Linear regression
    % Linear trends of 6 motion parameters
    % First order drifts of 6 motion parameters
    % Mean signal intensity of voxels in WM and CSF per segmentation
    
    motionParameters = dlmread(motionParametersFile)'; % TODO correct orientation
    nMotionParameters = size(motionParameters, 1);
    
    regressors = zeros(2*nMotionParameters + nMaskRegressors + GSRFlag, nTimePoints);
    
    % Linear trends and first order drifts of ~6 motion parameters.
    regressors(1:nMotionParameters, :) = motionParameters;
    regressors(nMotionParameters+1:2*nMotionParameters, 2:end) = ...
        diff(motionParameters, 1, 2); % first column are zeros
    
    % Mean signal intensity in regressionMask regions.
    for iR = 1:nMaskRegressors
        regressionIndx = segmentation == regressionMask(iR);
        regressors(2*nMotionParameters+iR, :) = ...
            mean(signalIntensities(regressionIndx(brainPrevalenceMask), :), 1);
    end
    
    % Global mean correction
    if GSRFlag
        regressors(end, :) = mean(signalIntensities, 1);
    end
    
    % Normalize regressors to avoid non-invertible matrices.
    regressors = bsxfun(@rdivide, regressors, mean(abs(regressors), 2));

    % Slow implementation:
    %     % Add column of constants for regress function.
    %     regressors = [ones(1, size(regressors, 2)); regressors];
    
    %     % Apply regressors
    %     selectedTimeSeries = fmri.Data.signalIntensities(selectedVoxels, :);
    %     
    %     for iV = 1:nnz(selectedVoxels)
    %         [x, xx, thisResiduals] = regress(signalIntensities(iV, :)', regressors');
    %         signalIntensities2(iV, :) = thisResiduals';
    %     end    
    
    % Fast implementation:
    % TODO: make this a function with appropriate tests.
    X = regressors';
    [Q,R,perm] = qr(X,0);
    X = X(:, perm);
    
    rc = rcond(R);
    if isnan(rc) || rc < 1e-6
        error(['Regression step unsuccessful. Covariate matrix is singular, ', ...
            'close to singular or badly scaled. RCOND = %g'], rc);
    end
    
    b = R \ (Q'* signalIntensities');
    residuals = signalIntensities - (X*b)'; % yhat = X*b
    
    signalIntensities = residuals - mean(residuals, 2);
    clear Q R perm X b residuals
    
    if pipelineV25Flag
        signalIntensities = signalIntensities(prevalenceMask, :); %#ok
        brainPrevalenceMask = false(size(brainMask));
        brainPrevalenceMask(brainMask) = prevalenceMask;
    end
    
    %% Bandpass filter
    
    % Get repetition time.
    props = getNiftiProperties(fmriProcessedFile);
    repetitionTimeMsec = props.res(4);
    assert(repetitionTimeMsec >= minRepetitionTime, ...
        ['Repetition time (%g msec) reported in fmriProcessedFile', ...
        ' is smaller than the minRepetitionTime (%g msec). ', ...
        'Transform repetition time to milliseconds ', ...
        'or adjust minRepetitionTime-parameter'], repetitionTimeMsec, minRepetitionTime);
    repetitionTimeSec = repetitionTimeMsec/1000;
    
    [filter_b, filter_a] = butter(2, 2*repetitionTimeSec*bandpass_filter.frequencies);
    
    % Use for-loop to avoid memory issues from having double. (1sec difference)
    % filteredSignal = filtfilt(filter_b, filter_a, double(signalIntensities(selectedVoxels, :)));
    filteredSignal = zeros(size(signalIntensities), 'single');
    for i = 1:size(signalIntensities, 1)
            filteredSignal(i,:) = filtfilt(filter_b, filter_a, ...
                double(signalIntensities(i, :)));
    end
    
    signalIntensities = filteredSignal;
    
    clear filteredSignal
    %% Scrub
    % Scrubbing removes frames from the rs-fMRI time-series that
    % display significant motion artifacts {Power, 2012 #50} before
    % correlation analysis.
    
    data = load(motionMetricsFile);
    motionMetrics = data.motionMetrics;
    metricDescriptions = data.metricDescriptions;
    clear data
    
    FD = motionMetrics(strcmp(metricDescriptions, 'FD'), :);
    DVARS = motionMetrics(strcmp(metricDescriptions, 'DVARS'), :);
    
    if scrubbing.scrubbing
        
        % When scrubbing is enabled (determined by the configuration
        % parameter scrubbing), frames with motion artifacts are identified
        % based on two indicators:
        
        % i) having framewise displacement FD larger than maxFD
        violations(1, :) = FD > scrubbing.maxFD;
        
        % ii) having a DVARS larger than Q3 + maxDVARS × IQR, where IQR
        % refers to the the interquartile range IQR = Q3 – Q1, with Q1 and
        % Q3 referring to the first and third quartile of the DVARS of all
        % frames.
        sortedDVARS = sort(DVARS, 'ascend');
        Q1 = sortedDVARS(round(0.25*length(sortedDVARS)));
        Q3 = sortedDVARS(round(0.75*length(sortedDVARS)));
        IQR = Q3 - Q1;
        violations(2, :) = DVARS > (Q3 + scrubbing.maxDVARS * IQR);
        
        % Frames with a number of indicators larger or equal to
        % minViolations are labeled as frames with potential motion
        % artifacts and are excluded from further analysis.
        outlierFrames = sum(violations, 1) >= scrubbing.minViolations;
        
        % To accommodate temporal smoothing of data, frames consecutive to
        % frames with labeled motion artifacts are optionally excluded:
        % configuration parameter backwardNeighbors determines the number
        % of preceding frames and forwardNeighbors determines the number of
        % succeeding frames to be excluded from further analysis.
        
        % This is a convolution, but for clarity let's use a for-loop
        outlierFramesextended = zeros(size(outlierFrames));
        for i = 1:length(outlierFrames)
            if outlierFrames(i) > 0
                outlierFramesextended(i-scrubbing.backwardNeighbors:i+scrubbing.forwardNeighbors) = 1;
            end
        end
        outlierFrames = outlierFramesextended;
    else
        outlierFrames = zeros(1, nTimePoints);
    end
    
    numberOfScrubbedVolumes = nnz(outlierFrames); %#ok
    
    % calculate average motion metric over remaining frames
    selectedFrames = find(~outlierFrames);
    motionMetrics = mean(motionMetrics(:, selectedFrames(2:end)), 2); %#ok
   
    % Compose final filtered, regressed and scrubbed time series
    selectedTimeSeries = signalIntensities(:, ~outlierFrames);
    clear signalIntensities
    
    %% Correlate   
    
    for iTemplate = 1:length(templates)
        thisTemplate = templates{iTemplate};
        fprintf('template: %s\n', thisTemplate);
        
        thisParcellationFile = strrep(parcellationFile, ...
            'TEMPLATE', thisTemplate);
        thisROIsFile = strrep(ROIsFile, ...
            'TEMPLATE', thisTemplate);
        thisConnectivityMatrixFile = strrep(strrep(connectivityMatrixFile, ...
            'METHOD', thisMethodDescription), ...
            'TEMPLATE', thisTemplate);
        thisTimeSeriesFile = strrep(strrep(timeSeriesFile, ...
            'METHOD', thisMethodDescription), ...
            'TEMPLATE', thisTemplate);
        thisLutFile = fullfile(templatesDir, thisTemplate, ...
            [thisTemplate '.annot.ctab']);
        
        parcellation = load_nifti(thisParcellationFile);
        parcellation = single(parcellation.vol(brainPrevalenceMask));
        ROIs = dlmread(thisROIsFile);
        
        % get associated regionDescriptions
        LUT = readtable(thisLutFile,  'filetype', 'text', ...
            'ReadVariableNames', false);
        LUT.Properties.VariableNames = {'ROIs', 'regionDescriptions', ...
            'Color1', 'Color2', 'Color3', 'Other'};
        [~, indxROIs] = ismember(ROIs, LUT.ROIs);
        regionDescriptions = deblank(LUT.regionDescriptions(indxROIs)); %#ok
        
        % Calculate time series averaged over all voxels of a region
        averageTimeSeries = zeros(length(ROIs), size(selectedTimeSeries, 2));
        for iR = 1:length(ROIs)
            averageTimeSeries(iR, :) = mean(selectedTimeSeries(...
                parcellation == ROIs(iR), :), 1);
        end
        
        if saveTimeSeriesFlag
            save(thisTimeSeriesFile, ...
                'averageTimeSeries', 'ROIs', 'regionDescriptions');
        end
        
        % Calculate correlation data
        [connectivity, pValues] = corr(averageTimeSeries');
        connectivity = connectivity .* ~eye(length(connectivity));
        pValues = pValues .* ~eye(length(connectivity)); %#ok
                
        save(thisConnectivityMatrixFile, ...
            'motionMetrics', 'numberOfScrubbedVolumes', ...
            'metricDescriptions', 'ROIs', ...
            'regionDescriptions', 'connectivity', 'pValues');
        
    end
    
    clear selectedTimeSeries
        
end

%% Clean up
status.reconstruction_functional_network = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_functional_network finished----\n');

