function [diffusionPeaks, GFA, QA] = gqi(signalIntensities, gtab, ...
    meanDiffusionDistanceRatio, outputPeaks, minPeakRatio, maxPeaks, nonlinearities)
% GQI     Reconstruct diffusion peaks through Generalized Q-ball Imaging
%
%   INPUT VARIABLES
%   signal_intensities:
%   2D matrix containing measured signal intensities.
%
%   gtab:
%   numberOfScansx3 matrix with the applied diffusion gradients as rows.
%   The norm of a gradient vector should be equal to the associated b-value.
%
%   meanDiffusionDistanceRatio:
%   Parameter regulating the coarseness of the reconstructed peak profile.
%   High values provide, theoretically, a more accurate reconstruction, but
%   also increase sensitivity to noise.
%
%   outputPeaks:
%   Maximum number of peaks per voxel included in the diffusion peaks file.
%
%   minPeakRatio:
%   Parameter controlling the sensitivity to detect peaks. Diffusion peaks
%   with a normalized coefficient (i.e. the coefficient of the peak divided
%   by the maximum coefficient) smaller than minPeakRatio are discarded
%
%   maxPeaks:
%   Number of identified peaks beyond which a voxel is considered
%   isotropic.
%
%   nonlinearities:
%   2D matrix containing the nonlinearity corrections per voxel.
%   Rows correspond to voxels and columns to corrections.
%
%   OUTPUT VARIABLES
%   diffusionPeaks:
%   Diffusion peak directions of each voxel. Diffusion
%   peaks are a Nx3xM matrix containing for N voxels at most M
%   diffusion peaks for each voxel. The first index corresponds to the
%   linear index of the voxel and the third index reflects the prominence
%   of the diffusion peak (the strongest peak having the lowest index). The
%   second dimension describes the direction of the diffusion peaks.
%
%   GFA:
%   Generalized fractional anisotropy measure of each voxel.
%
%   QA:
%   Quantitative anisotropy measure of each voxel.
%
%   NOTES
%   This function (following DSI Studio), uses 0.015 for 6D and thus the
%   diffusion length ratio appears longer. (e.g., a ratio 1.25 is in fact
%   1.145).

%   This code is adapted from the GQI implementation written by Fang-Cheng
%   (Frank) Yeh presented on:
%       http://dsi-studio.labsolver.org/Manual
%       /Reconstruction#TOC-Generalized-Q-sampling-Imaging-GQI-
%   gqi_reco.m:
%       http://docs.google.com/viewer?a=v&pid=sites&srcid=
%       bGFic29sdmVyLm9yZ3xkc2ktc3R1ZGlvfGd4OjE0YjBlMDQ5ODk5ODcwNTA
%   find_peak.m:
%       https://docs.google.com/viewer?a=v&pid=sites&srcid=
%       bGFic29sdmVyLm9yZ3xkc2ktc3R1ZGlvfGd4OjFmMGI4MzFmMWE0M2FmNjQ


if nargin < 7
    nonlinearitiesFlag = false;
else
    nonlinearitiesFlag = true;
end

dataDir = fullfile(fileparts(mfilename('fullpath')), 'reconstruction_basis');
reconBasis = dlmread(fullfile(dataDir, 'reconstruction_basis.txt'));
nPointsBasis = size(reconBasis, 1)/2;
reconBasis = reconBasis(1:nPointsBasis, :);

faces = dlmread(fullfile(dataDir, 'faces.txt'));
faces = faces - (faces > nPointsBasis)*nPointsBasis;

nVoxels = size(signalIntensities, 1);

if nonlinearitiesFlag
        
    G = sqrt(gtab.bvals) .* gtab.bvecs;
    odf = zeros(nVoxels, nPointsBasis);

    for iV = 1:nVoxels    
        bVec = sqrt(0.01506) .* G * (eye(3) + reshape(nonlinearities(iV, :), [3 3]));
        basisFunctions = meanDiffusionDistanceRatio .* bVec * reconBasis';
        basisFunctions = sin(basisFunctions) ./ basisFunctions;
        basisFunctions(isinf(basisFunctions) | isnan(basisFunctions)) = 1;
        
        odf(iV, :) = signalIntensities(iV, :) * basisFunctions;
    end
    
else
    
    bVec = sqrt(gtab.bvals*0.01506) .* gtab.bvecs;
    basisFunctions = meanDiffusionDistanceRatio .* bVec * reconBasis';
    basisFunctions = sin(basisFunctions) ./ basisFunctions;
    basisFunctions(isinf(basisFunctions) | isnan(basisFunctions)) = 1;
    
    odf = signalIntensities * basisFunctions;
    
end

GFA = sqrt((nPointsBasis * sum((odf - mean(odf, 2)).^2, 2)) ./ ...
    ((nPointsBasis-1) * sum(odf.^2, 2)));

peakIndices = zeros(size(signalIntensities, 1), outputPeaks, 'uint16');
QA = zeros(size(signalIntensities, 1), outputPeaks);

for i = 1:nVoxels
    
    thisOdf = odf(i, :);
    odfPeaks = thisOdf;
    odfPeaks(faces(thisOdf(faces(:, 2)) >= thisOdf(faces(:, 1)) | ...
        thisOdf(faces(:, 3)) >= thisOdf(faces(:, 1)), 1)) = 0;
    odfPeaks(faces(thisOdf(faces(:, 1)) >= thisOdf(faces(:, 2)) | ...
        thisOdf(faces(:, 3)) >= thisOdf(faces(:, 2)), 2)) = 0;
    odfPeaks(faces(thisOdf(faces(:, 1)) >= thisOdf(faces(:, 3)) | ...
        thisOdf(faces(:, 2)) >= thisOdf(faces(:, 3)), 3)) = 0;
    
    detectedPeaks = ( odfPeaks ./ max(odfPeaks) ) > minPeakRatio;
    
    % skip isotropic voxels
    numberOfPeaks = sum(detectedPeaks);
    if numberOfPeaks > maxPeaks
        continue;
    end
    
    % convert peak values into peak indices
    odfPeaks = odfPeaks .* detectedPeaks;
    
    if any(odfPeaks)
        [~, I] = sort(odfPeaks, 'descend');
        
        m = min(outputPeaks, nnz(odfPeaks));
        peakIndices(i, 1:m) = I(1:m);
        
        % Find GA of the peaks QA at u = Z(ODF(u)-iso(ODF)) the minimum value
        % of an ODF function is used as the isotropic component Z is a scaling
        % constant that makes the maximum of all iso(ODF) in the space equal to
        % one.
        QA(i, 1:m) = (odfPeaks(I(1:m)) - min(thisOdf)) ./ max(thisOdf);
    end
    
end

% find peaks
reconBasis = [0 0 0; reconBasis];
diffusionPeaks = zeros(size(peakIndices, 1), 3, size(peakIndices, 2));
for iP = 1:size(peakIndices, 2)
    diffusionPeaks(:, :, iP) = reconBasis(peakIndices(:, iP)+1, :);
end



