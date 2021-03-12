function diffusionTensor = iRESTORE(signalIntensities, gtab, thresCondNum, thresVarProjScores, nonlinearities)
% IRESTORE   Estimate diffusion tensor using iRESTORE algorithm.
%
%   iRESTORE models the signal of a voxel by a tensor and estimates one
%   preferred diffusion-direction per voxel, while reducing the impact of
%   physiological noise artifacts on the DTI modeling.
%
%   INPUT VARIABLES
%   signalIntensities:
%   2D matrix containing measured signal intensities.
%   Rows correspond to voxels and columns to diffusion gradients.
%
%   gtab:
%   numberOfScansx3 matrix with the applied diffusion gradients as rows.
%   The norm of a gradient vector should be equal to the associated b-value.
%
%   thresCondNum:
%   Condition number threshold.
%
%   thresVarProjScores:
%   projection variation threshold.
%
%   nonlinearities:
%   2D matrix containing the nonlinearity corrections per voxel.
%   Rows correspond to voxels and columns to corrections.   
%
%   OUTPUT VARIABLES
%   diffusionTensor:
%   nVoxels x 6 matrix with the estimated diffusion tensor parameters.

%   References:
%
%   de Reus, M. A. (2015). An eccentric perspective on brain networks
%   (Doctoral dissertation, Uitgeverij BOXPress).
%
%   Kingsley, P. B. (2006). Introduction to diffusion tensor imaging
%   mathematics: Part III. Tensor calculation, noise, simulations, and
%   optimization. Concepts in Magnetic Resonance Part A, 28(2), 155-179.
%
%   Chang, L. C., Jones, D. K., & Pierpaoli, C. (2005). RESTORE: robust
%   estimation of tensors by outlier rejection. Magnetic Resonance in
%   Medicine: An Official Journal of the International Society for Magnetic
%   Resonance in Medicine, 53(5), 1088-1095.
%
%   Chang, L. C., Walker, L., & Pierpaoli, C. (2012). Informed RESTORE: a
%   method for robust estimation of diffusion tensor from low redundancy
%   datasets in the presence of physiological noise artifacts. Magnetic
%   resonance in medicine, 68(5), 1654-1663.



MIN_POSITIVE_SIGNAL = 0.0001;
MAX_ITERATION = 400;

if nargin < 5
    nonlinearitiesFlag = false;
else 
    nonlinearitiesFlag = true;
end

signalIntensities = double(signalIntensities); % double is faster than single

nVoxels = size(signalIntensities, 1);
weightedScansAll = gtab.bvals > 0;

% Get B matrix with average b0-scans.
[gtabCorrected, Ball] = getGradientMatrix(gtab);

diffusionTensor = nan(nVoxels, 6, 'single');

warningsOld = warning();
warning('error', 'MATLAB:nearlySingularMatrix'); %#ok
warning('error', 'MATLAB:singularMatrix'); %#ok
warning('error', 'MATLAB:illConditionedMatrix'); %#ok
exitErrors = categorical(repmat({''}, nVoxels, 1));

for iVoxel = 1:nVoxels
    
    Sall = signalIntensities(iVoxel, :)';
    indx_succesful = Sall >= MIN_POSITIVE_SIGNAL;
    
    Sweighted = Sall(weightedScansAll & indx_succesful);
    S0 = Sall(~weightedScansAll & indx_succesful);
    
    if isempty(S0) || isempty(Sweighted)
        exitErrors(iVoxel) = 'b0-scan or all weighted scans unsuccessful.';
        continue;
    end
    
    % Following de Reus (2015), the six diffusion tensor parameters are
    % iteratively fitted, but the S0 signal is estimated by the average
    % signal of the unweighted scans.
    S = Sweighted ./ mean(S0);
    
    if nonlinearitiesFlag
        [gtabCorrected, Ball] = getGradientMatrix(gtab, nonlinearities(iVoxel, :));
    end
    
    B = Ball(indx_succesful(weightedScansAll), :);
    nWeighted = size(B, 1);
    
    % First a rough estimation is made of the diffusion tensor alphaVec
    % using a weighted linear least squares fit of the natural logarithm of
    % the signal.
    %
    % SigmaInv = diag((S.^2) / SigmaSquared) where SigmaSquared is the
    % same in each image (assuming all images have same noise level) and
    % disappears (eq 52, Kingsley 2005).
    
    x = log(S);
    SigmaInv = diag(S.^2);
    
    try
        alphaVec = (B'*SigmaInv*B)\(B'*SigmaInv)*x; % eq[51]
    catch ME
        if any(strcmp(ME.identifier, {'MATLAB:nearlySingularMatrix', ...
                'MATLAB:singularMatrix', 'MATLAB:illConditionedMatrix'}))
            exitErrors(iVoxel) = ['Not enough information in scan data ', ...
                'to make initial linear LS fit.'];
            continue
        else
            rethrow(ME)
        end
    end
    
    % The diffusion tensor is fitted using a nonlinear least-squares method
    % using the linear fit as initial guess of parameters.
    [alphaVec, ~, r, exitFlag] = lmsolverDTI(alphaVec, B, S);
    if exitFlag == 4
        exitErrors(iVoxel) = ['Nonlinear least-squares method did not ', ...
            'converge to initial guess of parameters.'];
       continue 
    end
    
    % Evaluate results of first fitting against goodnes-of-fit criterion.
    % If the residuals of all datapoints are within the CI (3*SD) then
    % assume no outliers and accept results with no further processing
    % (p1089, Chang 2004). SD is defined in (p1656, Chang 2012).
    %
    % TODO FEATURE: select voxels from centrum semiovale rostral to the corpus
    % callosum
    % TODO FEATURE: exclude some voxels as outliers
    SD = sqrt(nWeighted / (nWeighted - 6)) * median(abs(r - median(r)));
    
    if all(abs(r) < 3*SD) || (SD < eps)
        diffusionTensor(iVoxel, :) = alphaVec;
        continue
    end
    
    % Initiale iterative reweighting using GMM weighting function
    % The weight (omega) for each data point is normalized to the average of all lthe
    % weighting factors to yield the maximum likelihood (p1089, Chang
    % 2004). The convergence criterionis  less than one percent change of
    % fitted tensor parameters.
    alphaVecOld = alphaVec;
    counter = 0;
    while (counter < MAX_ITERATION)
        
        C = 1.4826 * median(abs(r - median(r)));
        omega = 1 ./ (r.^2 + C.^2);
        omega = omega ./ mean(omega);  
        
        if any(isnan(omega))
            break;
        end
        
        [alphaVec, ~, r, exitFlag] = lmsolverDTI(alphaVec, B, S, omega);
        
        if exitFlag == 4
            break
        end
        
        if all(abs((alphaVec - alphaVecOld)./alphaVecOld) < 0.01)
            break;
        end
        
        counter = counter + 1;
        alphaVecOld = alphaVec;

    end

    % Points outside the confidence interval (3*SD original signal) are
    % outliers (p1089, Chang 2004).
    SD = sqrt(nWeighted / (nWeighted - 6)) * median(abs(r - median(r)));
    indxOutliers = abs(r) > 3*SD;
    
    % Check b-matrix ill-conditioned or directionally unbalanced.
    % Coefficient of variation in the average projection scores calculated
    % following (eq. [2.10], page 38, de Reus 2015).
    condNum = cond(B(~indxOutliers, :));
    
    G = gtabCorrected.bvecs(weightedScansAll & indx_succesful, :);
    projScores = mean(abs(gtabCorrected.bvecs * G(~indxOutliers, :)'), 2);
    varProjScores = std(projScores) / mean(projScores);
    
    if condNum >= thresCondNum || varProjScores >= thresVarProjScores
        indxOutliers = false(size(indxOutliers));
    end
    
    % Final tensor fit with equal weights
    [alphaVec, ~, ~, exitFlag] = lmsolverDTI(alphaVec, ...
        B(~indxOutliers, :), S(~indxOutliers), ...
        1/SD^2 .* ones(nnz(~indxOutliers), 1));
    
    if exitFlag == 4
        exitErrors(iVoxel) = ['Final nonlinear least-squares method ', ...
            'did not converge.'];
    end
    
    diffusionTensor(iVoxel, :) = alphaVec;
    
end

if any(~ismissing(exitErrors))
   fprintf('Warning:\n');
   exitMessages = categories(exitErrors);
   exitCounts = countcats(exitErrors);
   for i=1:length(exitMessages)
       fprintf('In %i voxels: %s\n', exitCounts(i), ...
           exitMessages{i});
   end
end

warning(warningsOld);

end

function [gtabC, B] = getGradientMatrix(gtab, nonlinearities)
% GETGRADIENTMATRIX convert gradient table to B matrix.
%
% Notes: nonlinearities is a 3x3 matrix.
% Notes: nonlinearities correction based on HCP appendix II:
% https://www.humanconnectome.org/storage/app/media/documentation/
%       data_release/October2012_Release_Appendix4.pdf

if nargin < 2
    nonlinearities = zeros(3);
else
    nonlinearities = reshape(nonlinearities, [3 3]);
end

% Code for voxel-wise correction of dMRI gradients
I = eye(3);
v = gtab.bvecs*(I+nonlinearities);
n = sqrt(sum(v.^2, 2));

gtabC.bvecs = v ./ n;
gtabC.bvecs(n == 0, :) = 0;
gtabC.bvals = n.^2.*gtab.bvals;

% B is inverse design matrix
% Design matrix or B matrix assuming Gaussian distributed tensor model
weightedScans = gtabC.bvals > 0;
B = nan(length(gtabC.bvals), 6);  % eq [2]
B(:, 1) = -gtabC.bvecs(:, 1) .* gtabC.bvecs(:, 1) .* 1 .* gtabC.bvals;   % Bxx
B(:, 2) = -gtabC.bvecs(:, 2) .* gtabC.bvecs(:, 2) .* 1 .* gtabC.bvals;   % Byy
B(:, 3) = -gtabC.bvecs(:, 3) .* gtabC.bvecs(:, 3) .* 1 .* gtabC.bvals;   % Bzz
B(:, 4) = -gtabC.bvecs(:, 1) .* gtabC.bvecs(:, 2) .* 2 .* gtabC.bvals;   % Bxy
B(:, 5) = -gtabC.bvecs(:, 1) .* gtabC.bvecs(:, 3) .* 2 .* gtabC.bvals;   % Bxz
B(:, 6) = -gtabC.bvecs(:, 2) .* gtabC.bvecs(:, 3) .* 2 .* gtabC.bvals;   % Byz

% Include only weighted scans as we divide by b0.
B = B(weightedScans, :);

end
