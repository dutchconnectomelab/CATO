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
%   nVoxels x 7 matrix with the estimated diffusion tensor parameters.

MIN_POSITIVE_SIGNAL = 0.0001;
MAX_ITERATION = 400;
BVAL_THRESHOLD = 10;

if nargin < 5
    nonlinearitiesFlag = false;
else 
    nonlinearitiesFlag = true;
end

% double-type signalIntensities is faster than single-type
signalIntensities = double(signalIntensities);

nVoxels = size(signalIntensities, 1);
weightedScansAll = gtab.bvals >= BVAL_THRESHOLD;

% All initial values are called All (later outliers are excluded)
% Get b-vecs and first scan is (average) b0.
Ball = getGradientMatrix(gtab, BVAL_THRESHOLD);
Sall = [mean(signalIntensities(:, ~weightedScansAll), 2), ...
    signalIntensities(:, weightedScansAll)]';

% Prepare iterating over all voxels.
diffusionTensor = nan(nVoxels, 7, 'single');

warningsOld = warning();
warning('error', 'MATLAB:nearlySingularMatrix'); %#ok
warning('error', 'MATLAB:singularMatrix'); %#ok
warning('error', 'MATLAB:illConditionedMatrix'); %#ok
exitErrors = categorical(repmat({''}, nVoxels, 1));

for iVoxel = 1:nVoxels
    
    indx_succesful = Sall(:, iVoxel) >= MIN_POSITIVE_SIGNAL;

    % If b0-scan is unsuccessful, skip this voxel.
    if ~indx_succesful(1)
        exitErrors(iVoxel) = 'b0-scan unsuccessful.';
        continue;
    end    
    
    S = Sall(indx_succesful, iVoxel);
    
    if nonlinearitiesFlag
        Ball = getGradientMatrix(gtab, BVAL_THRESHOLD, nonlinearities(iVoxel, :));
    end
    
    B = Ball(indx_succesful, :);
    weightedScans = weightedScansAll(indx_succesful);
    nWeighted = nnz(any(B(:, 1:6), 2));
    
    
    % 1. A rough estimation of the diffusion tensor D is obtained using a
    % weighted linear least squares fit of the natural logarithm of the signal.
    
    % Erors will be less if the tensor elements are combined before calculating
    % the tensor.
    x = log(S); % eq [18]
       
    % SigmaInv = diag((S.^2) / SigmaSquared)
    % SigmaSquared is the same in each image (assuming all images have same
    % noise level) and disappears.
    SigmaInv = diag(S.^2); % eq [52]
    
    % alpha is seven-element column vector reconstruction of signal intensity
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
    
    % 2. The diffusion tensor is fitted using a nonlinear least-squares method
    % using the linear fit as initial guess of parameters.
    [alphaVec, resNorm, r, exitFlag] = lmsolverDTI(alphaVec, B, S);
    if exitFlag == 4
        exitErrors(iVoxel) = ['Nonlinear least-squares method did not ', ...
            'converge to initial guess of parameters.'];
       continue 
    end
    
    % Evaluate results of first fitting against goodnes-of-fit criterion
    
    % Test whether all data points lie within the confidence interval of 3xSD
    % of signal.
    % RMAD from iRESTORE paper p. 1656
    % TODO FEATURE: select voxels from centrum semiovale rostral to the corpus
    % callosum
    % TODO FEATURE: exclude some voxels as outliers
    CI = sqrt(nWeighted / (nWeighted - 7)) * median(abs(r - median(r)));
    
    % If the residuals of all datapoints are within the CI then assume no
    % outliers and accept results with no further processing.
    if all(abs(r) < 3*CI) || (CI < eps)
        diffusionTensor(iVoxel, :) = alphaVec;
        continue
    end
    
    % Initiale iterative reweighting using GMM weighting function
    % The weight for each data point is normalized to the average of al lthe
    % weighting factors to yeald the maximum likelihood
    
    % 3. An iterative reweighting process of the signal variance and outlier
    % detection is performed providing the final parameter fit.
    resNormOld = resNorm;
    counter = 0;
    while (counter < MAX_ITERATION)
        
        C = 1.4826 * median(abs(r - median(r)));
        omega = 1 ./ (r.^2 + C.^2);
        omega = omega ./ mean(omega);  % DIPY: weights normalized to mean weight (see p. 1089)
        
        if any(isnan(omega))
            break;
        end
        
        [alphaVec, resNorm, r, exitFlag] = lmsolverDTI(alphaVec, B, S, omega);
        
        if exitFlag == 4
            break
        end
        
        % convergence criterion: less than one percent increase
        if (abs(resNorm - resNormOld) / resNormOld) < 0.01
            break;
        end
        
        counter = counter + 1;
        resNormOld = resNorm;
    end

    
    % Points outside CI are identified as outliers and excluded
    CI = sqrt(nWeighted / (nWeighted - 7)) * median(abs(r - median(r)));
    indxOutliers = abs(r) > 3*CI;
    
    % sigmaSquared is the signal SD.
    % m = size(S, 1);
    % sigmaSquared = (s.^2) / (m - 6 + 1);
    
    % Condition number of matrix H (B) is smaller than threshold t_c after
    % eliminating outliers.
    condNum = cond(B(~indxOutliers, 1:6));
    
    % Calculate coefficient of variation in the average projection scores
    % (eq. [2.10], page 38)
    Bnorm = B ./ sqrt(sum(B.^2, 2));
    projScores = mean(abs(Bnorm(weightedScans, :) * Bnorm(weightedScans & ~indxOutliers, :)'), 2);
    varProjScores = std(projScores) / mean(projScores);
    
    if condNum >= thresCondNum || varProjScores >= thresVarProjScores
        indxOutliers = false(size(indxOutliers));
    end
    
    % Other points are weighted equally in nonlinear LS method
    S = S(~indxOutliers);
    B = B(~indxOutliers, :);
    
    [alphaVec, ~, ~, exitFlag] = lmsolverDTI(alphaVec, B, S);
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

function B = getGradientMatrix(gtab, bval_threshold, nonlinearities)
% GETGRADIENTMATRIX convert gradient table to B matrix.
%
% Notes: nonlinearities is a 3x3 matrix.
% Notes: nonlinearities correction based on HCP appendix II:
% https://www.humanconnectome.org/storage/app/media/documentation/
%       data_release/October2012_Release_Appendix4.pdf

if nargin < 3
    nonlinearities = zeros(3);
else
    nonlinearities = reshape(nonlinearities, [3 3]);
end

% Code for voxel-wise correction of dMRI gradients
I = eye(3);
v = gtab.bvecs*(I+nonlinearities);
n = sqrt(sum(v.^2, 2));

gtab.bvecs = v ./ n;
gtab.bvecs(n == 0, :) = 0;
gtab.bvals = n.^2.*gtab.bvals;

% B is inverse design matrix
% Design matrix or B matrix assuming Gaussian distributed tensor model
weightedScans = gtab.bvals > bval_threshold;
B = nan(length(gtab.bvals), 7);  % eq [2]
B(:, 1) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 1) .* 1 .* gtab.bvals;   % Bxx
B(:, 2) = -gtab.bvecs(:, 2) .* gtab.bvecs(:, 2) .* 1 .* gtab.bvals;   % Byy
B(:, 3) = -gtab.bvecs(:, 3) .* gtab.bvecs(:, 3) .* 1 .* gtab.bvals;   % Bzz
B(:, 4) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 2) .* 2 .* gtab.bvals;   % Bxy
B(:, 5) = -gtab.bvecs(:, 1) .* gtab.bvecs(:, 3) .* 2 .* gtab.bvals;   % Bxz
B(:, 6) = -gtab.bvecs(:, 2) .* gtab.bvecs(:, 3) .* 2 .* gtab.bvals;   % Byz
B(:, 7) = ones(size(gtab.bvals));

% First scan is (average) b0-weighted scan
B0 = mean(B(~weightedScans, :), 1);
Bweighted = B(weightedScans, :);
B = [B0; Bweighted];

end
