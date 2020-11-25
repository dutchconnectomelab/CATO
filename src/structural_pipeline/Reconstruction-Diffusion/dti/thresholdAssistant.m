function [thresCondNum, thresVarProjScores] = thresholdAssistant(gtab)
% THRESHOLDASSISTANT    Estimate iRESTORE thresholds using bootstrapping.
%
%   INPUT VARIABLES
%   gtab:
%   numberOfScansx3 matrix with the applied diffusion gradients as rows.
%   The norm of a gradient vector should be equal to the associated b-value.
%
%   OUTPUT VARIABLES
%   thresCondNum:
%   Suggested threshold for the condition number of the B-matrix.
%
%   thresVarProjScores:
%   Suggested threshold on the variation in the average projection scores.
%
%   NOTES
%   The condition number and variation in average projection scores
%   thresholds are specific for each gradient acquisition scheme. This
%   function estimates both thresholds using a bootstrapped sample of
%   gradient schemes obtained by randomly removing gradient directions from
%   the original scheme. Following experience from practice, as described
%   in de Reus (2015), thresholds where chosen from the distribution such
%   that:
%       - The removal of 5% of the gradients was allowed for 75% of the
%       samples.
%       - And the removal of 50% of the gradients was allowed in 25% of the
%       samples.
%   For each threshold these two points were obtained across all
%   permutations and the final thresholds were obtained by averaging the
%   two points.

%   Based on:
%   de Reus, M. A. (2015). An eccentric perspective on brain networks,
%   Uitgeverij BOXPress.
%   Chang, L. C., et al. (2012). "Informed RESTORE: A method for robust
%   estimation of diffusion tensor from low redundancy datasets in the
%   presence of physiological noise artifacts." Magnetic Resonance in
%   Medicine 68(5): 1654-1663.

%% Initialization

% Define B and G.
weightedScans = gtab.bvals > 0;

B = nan(length(gtab.bvals), 7);  % eq [2]
B(:, 1) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 1) .* 1 .* gtab.bvals;   % Bxx
B(:, 2) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 2) .* 2 .* gtab.bvals;   % Bxy
B(:, 3) = -gtab.bvecs(:, 2) .* gtab.bvecs(:, 2) .* 1 .* gtab.bvals;   % Byy
B(:, 4) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 3) .* 2 .* gtab.bvals;   % Bxz
B(:, 5) = -gtab.bvecs(:, 2) .* gtab.bvecs(:, 3) .* 2 .* gtab.bvals;   % Byz
B(:, 6) = -gtab.bvecs(:, 3) .* gtab.bvecs(:, 3) .* 1 .* gtab.bvals;   % Bzz
B(:, 7) = ones(size(gtab.bvals));

B = B(weightedScans, :);

G = gtab.bvecs;
G = G(weightedScans, :);

%% Bootstrap samples with randomly removed directions

nDir = nnz(weightedScans);
nPerm = 100;

nDirRemoved = nDir * [5 50] / 100;
nDirRemoved = max(round(nDirRemoved), 1);
nRemoved = length(nDirRemoved);

thresPercentages = [75 25];

condNum = nan(nPerm, nRemoved);
varProjScores = nan(nPerm, nRemoved);
for iPerm = 1:nPerm
    for iRemoved = 1:nRemoved
        
        indxDirRemain = randperm(nDir, nDir-nDirRemoved(iRemoved));
        
        condNum(iPerm, iRemoved) = cond(B(indxDirRemain, 1:6));
        
        projScores = mean(abs(G * G(indxDirRemain, :)'), 2);
        varProjScores(iPerm, iRemoved) = std(projScores) / mean(projScores);
        
    end
end

%% Select threshold based on random sample
thresCondNum = nan(nRemoved, 1);
thresVarProjScores = nan(nRemoved, 1);
for iR = 1:nRemoved
    
    thisThreshold = round(nPerm * thresPercentages(iR) / 100);
    
    condNum(:, iR) = sort(condNum(:, iR), 'descend');
    thresCondNum(iR) = condNum(thisThreshold, iR);
    
    varProjScores(:, iR) = sort(varProjScores(:, iR), 'descend');
    thresVarProjScores(iR) = varProjScores(thisThreshold, iR);
    
end

thresCondNum = mean(thresCondNum);
thresVarProjScores = mean(thresVarProjScores);