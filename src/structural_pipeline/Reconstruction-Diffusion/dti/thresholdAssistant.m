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
%   To ensure that enough information is preserved for reliable tensor
%   estimation, the iRESTORE algorithm (Chang, 2012) checks that the
%   B-matrix remains well conditioned and directionally balanced. The
%   B-matrix is considered well-conditioned if the condition number is lower
%   than thresCondNum and directionally balanced if the variation in average
%   projection scores is lower than thresVarProjScores. The condition number
%   and variation in average projection scores thresholds are specific for
%   each gradient acquisition scheme.
%
%   This function estimates each threshold as the average of two threshold
%   estimates (de Reus, 2015). Each threshold estimate is obtained from
%   bootstrapping a sample of gradient schemes in which random gradient
%   directions are removed from the original scheme (1,000 permutations).
%   
%      - The first threshold estimate is the value such that the removal of
%        5% of the gradients is accepted in 75% of the samples.
%      - The second threshold estimate is the value such that the removal
%        of 50% of the gradients is  accepted in 25% of the samples.
%
%   The suggested threshold for the condition number and variation in
%   average projection scores is then the average of these two estimates.

%   Based on:
%   de Reus, M. A. (2015). An eccentric perspective on brain networks,
%   Uitgeverij BOXPress.
%   Chang, L. C., et al. (2012). "Informed RESTORE: A method for robust
%   estimation of diffusion tensor from low redundancy datasets in the
%   presence of physiological noise artifacts." Magnetic Resonance in
%   Medicine 68(5): 1654-1663.

%% Initialization

NREMOVED = [5 50];
PERCENTAGE_ALLOWED = [75 25];
NPERM = 1000;

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
nDirRemoved = nDir * NREMOVED / 100;
nDirRemoved = max(ceil(nDirRemoved), 1);

nPoints = length(NREMOVED);
thresCondNum = nan(nPoints, 1);
thresVarProjScores = nan(nPoints, 1);

for iRemoved = 1:nPoints
    
    condNum = nan(NPERM, 1);
    varProjScores = nan(NPERM, 1);
    
    for iPerm = 1:NPERM
        indxDirRemain = randperm(nDir, nDir-nDirRemoved(iRemoved));
        
        condNum(iPerm) = cond(B(indxDirRemain, 1:6));
        
        projScores = mean(abs(G * G(indxDirRemain, :)'), 2);
        varProjScores(iPerm) = std(projScores) / mean(projScores); 
    end
    
    thisThreshold = round(NPERM * PERCENTAGE_ALLOWED(iRemoved) / 100);
    
    condNum = sort(condNum, 'ascend');
    thresCondNum(iRemoved) = condNum(thisThreshold);
    
    varProjScores = sort(varProjScores, 'ascend');
    thresVarProjScores(iRemoved) = varProjScores(thisThreshold);    
    
end

thresCondNum = mean(thresCondNum);
thresVarProjScores = mean(thresVarProjScores);