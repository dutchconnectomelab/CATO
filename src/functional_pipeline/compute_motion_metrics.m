function compute_motion_metrics(configParams)

%% Initialization
status.compute_motion_metrics = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---Compute_motion_metrics started----\n');

%% Computation motion metrics
motionParametersFile = configParams.functional_preprocessing.motionParametersFile;
fmriProcessedFile = configParams.functional_preprocessing.fmriProcessedFile;
motionMetricsFile = configParams.compute_motion_metrics.motionMetricsFile;
segmentationFile = configParams.functional_preprocessing.segmentationFile;

motionParams = dlmread(motionParametersFile);

props = getNiftiProperties(fmriProcessedFile);
ntimepoints = props.dim(4);

metricDescriptions = {'FD'; 'DVARS'};
motionMetrics = zeros(length(metricDescriptions), ntimepoints);

% compute framewise displacement (FD)
% framewise displacement (FD) is defined as the sum of the estimated
% translational and rotational displacement in a frame. Rotational
% displacement is modeled on a sphere of radius 50 mm
radius = 50; % mm
motionMetrics(1, 2:end) = ...
    abs(motionParams(2:end, :) - motionParams(1:end-1,:)) * ...
    [radius radius radius 1 1 1]';

% DVARS is calculated as the square root of the average squared
% intensity differences of all brain voxels between two consecutive
% frames in the rs-fMRI volume
% Create brain mask
segmentation = load_nifti(segmentationFile);
maskBrain = segmentation.vol(:) > 0;

signalIntensities = load_nifti_fmri(configParams, maskBrain);
maskSignal = any(signalIntensities,2);
signalMaskedDifference = signalIntensities(maskSignal, 2:end) - ...
    signalIntensities(maskSignal, 1:end-1);

clear signalIntensities

motionMetrics(2, 2:end) = sqrt(mean(signalMaskedDifference.^2, 1)); %#ok

save(motionMetricsFile, 'metricDescriptions', 'motionMetrics');

%% Clean up
status.compute_motion_metrics = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---Compute_motion_metrics finished----\n');

