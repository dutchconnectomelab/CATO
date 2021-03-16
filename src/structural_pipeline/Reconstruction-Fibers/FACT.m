function fibers = FACT(diffusionPeaks, seedLocations, mask, voxelSize, maxAngleDeg, maxFiberRadius)
% FACT Reconstruct fibers using Fiber Assignment by Continuous Tracking
%
% INPUT VARIABLES
%   diffusionPeaks:
%   Diffusion peak directions of each voxel. Diffusion peaks are a Nx3xM
%   matrix containing for N voxels at most M diffusion peaks for each
%   voxel.
%
%   mask:
%   Voxels in the mask with values 0 are not entered, -1 are stopping
%   points and all other voxels have value 1.
%
%   voxelSize:
%   Voxel size of dwiProcessedFile.
%
%   maxAngleDeg:
%   Largest turn in degrees a fiber is allowed to take. Fiber
%   reconstruction stops if a tracker is about to make a sharp turn (with
%   angle > maxAngleDeg).
%
%   maxFiberRadius:
%   Maximum number of steps from seed to endpoints. Fiber reconstruction
%   stops if the number of steps from the seed is larger than
%   maxFiberRadius. (Maximum length of fibers in mm depends on the voxel
%   size).
%
%   OUTPUT VARIABLES
%   fibers:
%   nSteps x 3 x nFibers tensor describing the reconstructed fibers.
%
%   NOTES
%   This function implements an extension of the “Fiber Assignment by
%   Continuous Tracking” (FACT) algorithm by allowing multiple peaks per
%   voxel.

%   REFERENCES
%   Mori, S., et al. (1999). "Three-dimensional tracking of axonal
%   projections in the brain by magnetic resonance imaging." Annals of
%   neurology 45(2): 265-269.

% INITIALIZE
minAngleArccos = cos(deg2rad(maxAngleDeg));
invVoxelSize = 1 ./ voxelSize;
nSeedLocations = size(seedLocations, 1);

assert(all(isreal(diffusionPeaks)), ['Diffusion peaks file ', ...
    'contains non-real entries. Make sure diffusion reconstruction step ', ...
    'completed sucessfully.']);

% mask:
% 3D matrix representing voxel space.
% Trackers will not enter voxels labeled by 0 and terminate if a voxel
% labeled by -1 is reached. Other voxels should be labeled by 1.

n1 = size(mask,1);
n2 = size(mask, 2);
n3 = size(mask, 3);

% Start trackers in seed voxels
% Trackers live in the continuous voxel space (writeFibers converts them to
% VOXMM space).
tracker = seedLocations; % continuous voxel space
trackerInd = floor(tracker);
trackerVoxelInd = trackerInd(:, 1)+(trackerInd(:,2)-1)*n1+(trackerInd(:,3)-1)*n1*n2;     
trackerVec = diffusionPeaks(trackerVoxelInd, :, 1); % VOXMM space.

% First half trackers goes forward, second half goes backwards (in parallel) 
tracker = [tracker; tracker];
trackerVec = [-trackerVec; trackerVec];

trackerIndx = 1:2*nSeedLocations;
previousMask = ones(2*nSeedLocations,1);

fibers = zeros(3, 4*maxFiberRadius+1, nSeedLocations, 'single');
fibers(:, 2*maxFiberRadius + 1, :) = seedLocations'; % start point of each tracker
fiberStepCount = zeros(2*nSeedLocations, 1); % stepwise length of each fiber
trackerIndPrevious = zeros(2*nSeedLocations, 2);
trackerIndPrevious(:, 1) = [trackerVoxelInd; trackerVoxelInd];


for iStep = 1:maxFiberRadius
        
    stopFlag = false(length(trackerIndx), 1);
    
    % Propagate trackers
    trackerInd = floor(tracker);
    trackerVoxelInd = trackerInd(:, 1)+(trackerInd(:,2)-1)*n1+(trackerInd(:,3)-1)*n1*n2; 
    
    trackerPeaks = diffusionPeaks(trackerVoxelInd, :, :);
    
    % look for best fit (trackerVec and diffusionPeaks are in VOXMM space).
    trackerPeakAngles = sum(trackerVec .* trackerPeaks, 2);
    [maxtrackerPeakAngle, I] = max(abs(trackerPeakAngles), [], 3);
    
    % Stop if turn > maxAngle
    stopFlag = stopFlag | (maxtrackerPeakAngle < minAngleArccos);
    
    % Get the peak number (first column, second column etc) for each voxel
    % of the peak that is closest to trackerVec. And update trackerVec
    trackerPeaks = permute(trackerPeaks, [1 3 2]);
    trackerPeaks = reshape(trackerPeaks, [], 3);
    minPeaksInd = [1:numel(trackerVoxelInd)]' + (I-1)*numel(trackerVoxelInd);
    trackerVec = trackerPeaks(minPeaksInd, :);
    
    % if angle is more than 90 degrees, flip the direction of the peak.
    n = size(trackerPeakAngles, 1);
    selectedPeakAngles = trackerPeakAngles(([1:n]' + (I-1)*n));
    trackerVec = trackerVec .* sign(selectedPeakAngles);
    
    % find voxel boundary
    trackerVecScaled = trackerVec .* invVoxelSize'; % get voxel space.
    boundaries = floor(tracker) + (trackerVecScaled >= 0); % FLits
    t = (boundaries - tracker) ./ trackerVecScaled;
    t = sort([t (t+1./abs(trackerVecScaled))], 2);
    tavg = 0.5*(t(:,1) + t(:,2)); % take average method.
    trackerBoundary = tracker + t(:,1).*trackerVecScaled;
    tracker = tracker + tavg.*trackerVecScaled;
    
    % Stop if outside mask
    stopFlag = stopFlag | any((floor(tracker) < 1) | (floor(tracker) > [n1 n2 n3]), 2);
    
    % Stop if this is a forbidden region (mask == 0)
    trackerInd = max(min(floor(tracker), [n1 n2 n3]), 1);
    trackerInd = trackerInd(:, 1)+(trackerInd(:,2)-1)*n1+(trackerInd(:,3)-1)*n1*n2;
    
    stopFlag = stopFlag | (mask(trackerInd) == 0);
    
    % Stop if previous region was a stop region (mask == -1)
    stopFlag = stopFlag | previousMask == -1;
    previousMask = mask(trackerInd);
    
    % Stop if tracker enters voxel that was visited in last 2 steps.
    stopFlag = stopFlag | any(trackerInd == trackerIndPrevious, 2);
    trackerIndPrevious = [trackerInd trackerIndPrevious(:,1)];
    
    % Stop tracking
    trackerIndx = trackerIndx(~stopFlag);
    tracker = tracker(~stopFlag, :);   
    trackerBoundary = trackerBoundary(~stopFlag, :);
    trackerVec = trackerVec(~stopFlag, :);
    previousMask = previousMask(~stopFlag);
    trackerIndPrevious = trackerIndPrevious(~stopFlag, :);
    
    % no trackers left.
    if isempty(trackerIndx)
        break
    end
    
    % Update fiber list
    % First set goes backwards second set goes forwards
    % First back is maxFiberRadius, first forward is maxFiberRadius+1
    indx = trackerIndx <= nSeedLocations;
    if any(indx)
        fibers(:, 2*maxFiberRadius - 2*iStep + 2, trackerIndx(indx)) = trackerBoundary(indx,:)';
        fibers(:, 2*maxFiberRadius - 2*iStep + 1, trackerIndx(indx)) = tracker(indx,:)';
        fiberStepCount(trackerIndx(indx)) = iStep;
    end
    if any(~indx) 
        fibers(:, 2*maxFiberRadius + 2*iStep, trackerIndx(~indx) - nSeedLocations) = trackerBoundary(~indx,:)';
        fibers(:, 2*maxFiberRadius + 2*iStep + 1, trackerIndx(~indx) - nSeedLocations) = tracker(~indx,:)';
        fiberStepCount(trackerIndx(~indx)) = iStep;
    end
end

% Exclude fibers traversing only 2 voxels.
% fiberStepCount = fiberStepCount(1:nSeedLocations) ...
%     + fiberStepCount(nSeedLocations+1:2*nSeedLocations);
% fibers = fibers(:, :, fiberStepCount>2);
