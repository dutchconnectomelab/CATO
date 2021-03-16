function [fiberProperties, propertyDescriptions] = ...
    getFiberPropertiesFromFile(fiberFile, ROIlist, parcellation, ...
    diffusionMeasures, weightDescriptions, includeGMVoxelsFlag)
% GETFIBERPROPERTIESFROMFILE    Calculate fiber properties of fiber file.
%
%   INPUT VARIABLES
%   fiberFile:
%   TRK file in which the reconstructed fiber cloud is saved. See:
%   http://trackvis.org/docs/?subsect=fileformat
%
%   ROIlist:
%   Parcellation codes of regions in connectivity matrix.
%
%   parcellation:
%   Structure (in NIFTI format) with brain parcellation.
%
%   diffusionMeasures:
%   nVoxels x nMeasures matrix with diffusion measures per voxel.
%
%   weightDescriptions:
%   Names of diffusion measures.
%
%   includeGMVoxelsFlag:
%   Flag to run fiber property reconstruction compatible with old CATO
%   versions.
%
%   OUTPUT VARIABLES
%   fiberProperties:
%   nSegments x nProperties matrix describing for each segment the fiber
%   properties.
%
%   propertyDescriptions:
%   List of properties in fiberProperties variable.

%% Initialization
if nargin < 5
    includeGMVoxelsFlag = false;
end

n1 = size(parcellation.vol,1);
n2 = size(parcellation.vol, 2);
nregions = length(ROIlist);

indexMatrix = zeros(nregions);
indexMatrix(tril(ones(nregions), -1) == 1) = 1:nregions*(nregions-1)/2;
indexMatrix = max(indexMatrix, indexMatrix');

maxEntries = 4e6; % based on an arbitrary limit, st. final matrix is 1.5GB
fiberIndex1 = single(zeros(maxEntries, 1));
fiberIndex2 = single(zeros(maxEntries, 1));
fiberMeasures = single(zeros(size(diffusionMeasures, 2), maxEntries));
fiberLength = single(zeros(maxEntries, 1));
fiberNumber = single(zeros(maxEntries, 1));
fiberMaxAngle = single(zeros(maxEntries, 1));
connectionNumber = single(zeros(maxEntries, 1));

% Prepare ROIlist and parcellationVol such that we can use ismembc (which
% is much faster than ismember).
ROIlist = unique(ROIlist);
ROIlist = single(ROIlist);
parcellationVol = single(parcellation.vol);

%% Loop through fibers
fid = fopen(fiberFile, 'r');
header = readTrkHeader(fid);
invVoxelSize = (1./ header.voxel_size)';

counter = 1;
iFiber = 0;
while true
    
    % Read fiber
    iFiber = iFiber + 1;
    points = fread(fid, 1, 'int');
    thisFiber = fread(fid, [(3 + header.n_scalars) points], 'float');
    
    if header.n_properties > 0
        fread(fid, header.n_properties, 'float');
    end
    
    if isempty(points)
        break;
    end
    
    % Skip fiber if there are less than 3 steps (2 GM voxels and 1 WM voxel).
    if points < 5
        continue
    end
    
    % Assert fibers are in point-boundary-point format.
    tmp = thisFiber(:, [2 end-1]) .* invVoxelSize;
    assert(all(any(abs(round(tmp) - tmp) < 1e-3, 1)), ...
        'TRK file must be in voxel-boundary format');
    
    % Find intersection with ROIs.
    fiberPoints = thisFiber(:, 1:2:end);
    fiberPointsInd = floor(fiberPoints .* invVoxelSize) + 1;
    fiberPointsInd = fiberPointsInd(1, :) + ...
                     (fiberPointsInd(2,:)-1)*n1 + ...
                     (fiberPointsInd(3,:)-1)*n1*n2;
    assert(min(fiberPointsInd(:)) > 0);
    assert(max(fiberPointsInd(:)) < numel(parcellationVol));
    fiberPointsROIs = parcellationVol(fiberPointsInd);
    
    % Skip fiber if no or only one ROI is touched.
    touchedROIs = ismembc(ROIlist, sort(fiberPointsROIs));
    if sum(touchedROIs) <= 1
        continue
    end
    
    % Calculate distance matrix with fiber distance between fiber points
    % in ROI regions.
    stepDistances = [0 sqrt(sum(diff(thisFiber,1,2).^2)) 0];
    stepDistancesPoints = stepDistances(1:2:end) + stepDistances(2:2:end);
    
    fiberPointsROIsInd = find(ismembc(fiberPointsROIs, ROIlist));
    roiPoints = fiberPointsROIs(fiberPointsROIsInd);
    
    distanceM = cumsum(stepDistancesPoints);
    distanceM = abs(distanceM(fiberPointsROIsInd) - distanceM(fiberPointsROIsInd)');
    
    % Calcualte fiber turn angles
    D = diff(fiberPoints,1,2);
    D = D ./ sum(D.^2,1).^(1/2);
    anglesNA = min(sum(D(:,2:end) .* D(:,1:end-1), 1), 1);
    anglesRadPoints = [0 acos(anglesNA) 0];
    
    % Calculate for each region pair that the fiber touches the fiber
    % segment properties.
    touchedROIsIndx = find(touchedROIs);
    touchedROIsCode = ROIlist(touchedROIs);
    for p = 1:length(touchedROIsIndx)
        for q = p+1:length(touchedROIsIndx)
            
            % Get fiber indices (where does it penetrate ROI q and ROI p).
            PIndx = find(roiPoints == touchedROIsCode(p));
            QIndx = find(roiPoints == touchedROIsCode(q));
            
            % Select shortest segment between the two ROIs.
            distanceMM = distanceM(PIndx, QIndx);
            [I, J] = find(distanceM(PIndx, QIndx) == min(distanceMM(:)),1);
            
            PtouchedIndx = fiberPointsROIsInd(PIndx(I));
            QtouchedIndx = fiberPointsROIsInd(QIndx(J));
            
            if PtouchedIndx > QtouchedIndx
                segmentIndx = QtouchedIndx+1:PtouchedIndx-1;
            elseif QtouchedIndx > PtouchedIndx
                segmentIndx = PtouchedIndx+1:QtouchedIndx-1;
            end
            
            % Compatibility option: exclude GM-GM fibers with no white
            % matter voxels.
            if ~includeGMVoxelsFlag
                if (abs(PtouchedIndx - QtouchedIndx) <= 1)
                    continue
                end
            end
            
            distancesSegmentIndx = stepDistancesPoints(segmentIndx);
            
            % Special case of GM - WM -GM with WM point on corner of voxel.
            if distancesSegmentIndx == 0
                continue
            end
            
            fiberLength(counter) = sum(distancesSegmentIndx);
            
            % Compatibility option.
            if includeGMVoxelsFlag
                if (abs(PtouchedIndx - QtouchedIndx) == 1)
                    % Include GM-GM fibers
                    segmentIndx = min(QtouchedIndx, PtouchedIndx): ...
                        max(QtouchedIndx, PtouchedIndx);
                    distancesSegmentIndx = ones(size(segmentIndx));
                    fiberLength(counter) = stepDistances(2*(segmentIndx(1)-1)+2) + ...
                        stepDistances(2*(segmentIndx(end)-1)+1);
                else
                    % Add GM-border length to fiber.
                    fiberLength(counter) = sum(distancesSegmentIndx) + ...
                        stepDistances(2*(segmentIndx(1)-1)) + ...
                        stepDistances(3+2*(segmentIndx(end)-1));
                end    
            end
            
            if counter > maxEntries
                error(['Number of reconstructed fiber segments exceeds ', ...
                    'maximum (maxEntries = %.g).'], ...
                    maxEntries);
            end
            
            % Prepare and save fiber segment measures.
            fiberMaxAngle(counter) = max(anglesRadPoints(segmentIndx));
            
            % Convert points to fiber indices (point-boundary-point)
            fiberIndex1(counter) = 1+2*(segmentIndx(1)-1);
            fiberIndex2(counter) = 1+2*(segmentIndx(end)-1);
            
            measuresSegment = diffusionMeasures(fiberPointsInd(segmentIndx), :);
            fiberMeasures(:, counter) = ...
                sum((distancesSegmentIndx' .* measuresSegment),1) ...
                ./ sum(distancesSegmentIndx);
            
            connectionNumber(counter) = indexMatrix(touchedROIsIndx(p), ...
                touchedROIsIndx(q));
            
            fiberNumber(counter) = iFiber;
            
            counter = counter + 1;
        end
    end
end

% Close reconstructed fiber file.
fclose(fid);

% Return fiberProperties and propertyDescriptions
fiberProperties = [fiberNumber(1:counter-1), ...
    connectionNumber(1:counter-1), ...
    fiberIndex1(1:counter-1), ...
    fiberIndex2(1:counter-1, :), ...
    fiberMaxAngle(1:counter-1), ...
    fiberLength(1:counter-1, :), ...
    fiberMeasures(:, 1:counter-1)'];
propertyDescriptions = {'fiberNumber', 'connectionNumber', ...
    'fiberIndex1', 'fiberIndex2', 'maxAngle' 'length', weightDescriptions{:}};
