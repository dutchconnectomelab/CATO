function reconstruction_fibers(configParams)
% RECONSTRUCTION_FIBERS   Reconstruct fibers from diffusion peaks.
%
%   reconstruction_fibers(CONFIGPARAMS) tracks fibers using the FACT
%   algorithm using the parameters as specified in CONFIGPARAMS.

%% Initialization

status.reconstruction_fibers = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---Reconstruction fibers started----\n');

% Prepare variables
maxMemoryGB = configParams.general.maxMemoryGB;
maxAngleDeg = configParams.reconstruction_fibers.maxAngleDeg;
maxFiberRadius = configParams.reconstruction_fibers.maxFiberRadius;
segmentationFile = configParams.structural_preprocessing.segmentationFile;
diffusionPeaksFile = configParams.reconstruction_diffusion.diffusionPeaksFile;
fiberFile = configParams.reconstruction_fibers.fiberFile;
diffusionMeasuresFile = configParams.reconstruction_diffusion.diffusionMeasuresFile;

forbiddenRegions = configParams.reconstruction_fibers.forbiddenRegions;
stopRegions = configParams.reconstruction_fibers.stopRegions;
startRegions = configParams.reconstruction_fibers.startRegions;
nSeedsPerVoxel = configParams.reconstruction_fibers.NumberOfSeedsPerVoxel;
minFA = configParams.reconstruction_fibers.minFA;
FAMeasure = configParams.general.FAMeasure;

reconMethods = configParams.general.reconstructionMethods;
if ~iscell(reconMethods)
    reconMethods = {reconMethods};
end

% Load segmentation
data = load_nifti(segmentationFile);
voxelSize = single(data.pixdim(2:4));
trkheader.voxel_size = data.pixdim(2:4);
trkheader.vox_to_ras = data.vox2ras;
segmentationVol = data.vol;
clear data

props = getNiftiProperties(segmentationFile);
orientation = props.orientation;

header = createTrkHeader(...
    'dim', size(segmentationVol), ...
    'voxel_size', trkheader.voxel_size, ...
    'vox_to_ras', trkheader.vox_to_ras, ...
    'voxel_order', orientation, ...
    'image_orientation_patient', [1 0 0 0 1 0]);

for iMethod = 1:length(reconMethods)
thisMethod = reconMethods{iMethod};

%% Initialize

disp([thisMethod ' fiber reconstruction started']);

% Prepare dynamic variables
thisDiffusionPeaksFile = strrep(diffusionPeaksFile, 'METHOD', thisMethod);
thisFiberFile = strrep(fiberFile, 'METHOD', thisMethod);
writeFibers([], thisFiberFile, voxelSize, header);

% Load dynamic data
data = load(thisDiffusionPeaksFile, 'diffusionPeaks');
diffusionPeaks = single(data.diffusionPeaks);
clear data

%% Prepare Mask
mask = zeros(size(segmentationVol)); % , 'int8'
mask(any(any(diffusionPeaks, 3), 2)) = 1;
mask(ismember(segmentationVol, stopRegions)) = -1;
mask(ismember(segmentationVol, forbiddenRegions)) = 0;

% Do not perform tracking in voxels with no diffusion peaks
indx = any(isnan(diffusionPeaks(:,:,1)),2);
mask(indx) = -1;

% Add voxels with FA < minFA as stop regions to mask.

% Dependent on the reconstruction method FA can be measured as FA or e.g.
% generalized fractional anisotropy.
FAMeasureMethods = fieldnames(FAMeasure);
indxFAMeasure = strcmpi(FAMeasureMethods, thisMethod);
if nnz(indxFAMeasure) == 0
    error('Cannot find reconstruction method (%s) in FAMeasure parameter (%s).', ...
        thisMethod, lower(jsonencode(FAMeasure)));
elseif nnz(indxFAMeasure) > 2
    error('Multiple method (%s) instances found in FAMeasure parameter (%s).', ...
        thisMethod, lower(jsonencode(FAMeasure)));
end

data = load(diffusionMeasuresFile);
indxFA = strcmpi(FAMeasure.(FAMeasureMethods{indxFAMeasure}), data.weightDescriptions);
assert(nnz(indxFA) == 1, ...
    'Cannot find fractional anisotropy measure (%s) in diffusionMeasuresFile (%s).', ...
    FAMeasure.(FAMeasureMethods{indxFAMeasure}), diffusionMeasuresFile);
maskFA = data.diffusionMeasures(:, :, :, indxFA) <= minFA;
mask(mask == 1 & maskFA) = -1; 
clear data

%% Prepare seed locations
% Get voxels with seeds
seedLocations = ismember(segmentationVol, startRegions) & (mask == 1);

% Get indices of seed voxels
[seedLocationI, seedLocationJ, seedLocationK] = ...
    ind2sub(size(seedLocations), find(seedLocations)); 
seedLocations = single([seedLocationI seedLocationJ seedLocationK]);

% Distribute seed locations evenly in voxels
m = nSeedsPerVoxel^(1/3);
assert(m == round(m));
[deltaI, deltaJ, deltaK] = meshgrid(1/(m+1):1/(m+1):1-1/(m+1));
deltaXYZ = [deltaI(:) deltaJ(:) deltaK(:)]; % in voxel space

% FEATURE: Use seeds per mm3?! Important for non-isotropic voxels
seedLocations = repmat(seedLocations, [1 1 nSeedsPerVoxel]) + permute(deltaXYZ, [3 2 1]);
seedLocations = permute(seedLocations, [3 1 2]);
seedLocations = reshape(seedLocations, [], 3);

%% Perform fiber tracking

    nSeeds = size(seedLocations, 1);
   
    % memory used in FACT algorithm:
    % (more-or-less) 3*(4*maxFiberRadius+1)*nSeedLocations*4
    % maxMemoryGB = maxMemoryGB*1e9 bytes
    nSeedsBatch = floor((maxMemoryGB*1e9) ./ (3*(4*maxFiberRadius+1)*4));
    
    for i = 1:nSeedsPerVoxel
        seedsInBatch = (i-1)*nSeedsBatch+1:min(i*nSeedsBatch, nSeeds);
        thisSeedLocations = seedLocations(seedsInBatch, :);
        
        fibers = FACT(diffusionPeaks, thisSeedLocations, mask, voxelSize, ...
            maxAngleDeg, maxFiberRadius);
        
        writeFibers(fibers, thisFiberFile, voxelSize)
        clear fibers

    end
    
disp([thisMethod ' fiber reconstruction finished']);

end

status.reconstruction_fibers = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---Reconstruction fibers finished----\n');

