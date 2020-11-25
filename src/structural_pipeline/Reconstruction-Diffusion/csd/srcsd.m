function [diffusion_peaks, test] = srcsd(signal_intensities, ...
    gtab, response, shOrder, lambda, tau, outputPeaks, ...
    minPeakRatio, maxPeaks, nonlinearities)
% SRCSD   Performs super resolved constrained spherical deconvolution.
%
%   These functions are based on the implementation of Dipy and described
%   in:
%     Tournier, J.D., et al. NeuroImage 2007. Robust determination of
%     the fibre orientation distribution in diffusion MRI:
%     Non-negativity constrained super-resolved spherical
%     deconvolution
% 
% INPUT VARIABLES
% signal_intensities:
% 2D matrix containing measured signal intensities.
% Rows correspond to voxels and columns to diffusion gradients.
%
% gradient_table:
% numberOfScansx3 matrix with the applied diffusion gradients as rows.
% The norm of a gradient vector should be equal to the associated b-value.
%
% response:
% Structure with 1) response.S0 describing the average intensity in b0
% scans and 2) response.evals a 3x1 vector with the first eigenvalue and
% average of second and third eigenvalue.
% 
% sh_order:
% Maximal spherical harmonics order.
%
% lambda:
% Weight given to the constrained-positivity regularization part of
% the deconvolution equation.
%
% tau:
% Threshold controlling the amplitude below which the corresponding
% fODF is assumed to be zero.  Ideally, tau should be set to
% zero. However, to improve the stability of the algorithm, tau is
% set to tau*100 % of the mean fODF amplitude.
%
% output_peaks:
% Maximum number of peaks per voxel in output.
%
% min_peak_ratio:
% Proportion of highest peak value that smaller peaks must exceed.
%
% max_peaks:
% Number of peaks beyond which a voxel is declared isotropic.
%
% OUTPUT VARIABLES
% peak_indices:
% 2D matrix listing for each voxel which vectors from the reconstruction 
% basis (see below) appear as local diffusion peaks.
% Rows correspond to voxels and column numbers indicate the relative size
% of the peaks (largest peaks are listed first).
%
% recon_basis_xyz:
% Mx3 matrix consisting of unit vectors, distributed uniformly over a half
% sphere, that may appear as diffusion peaks.
% Peak indices correspond to rows of this matrix.
%
% test:
% Variable with intermediate output for debugging and test prupose.

%%  test environment
% set test.status to true to get intermediate output for comparison with
% dipy implementation.

test.status = false; % no test statistics

%% Initializatie van iets anders

if nargin < 10
    nonlinearitiesFlag = false;
else
    nonlinearitiesFlag = true;
end

warningsOld = warning();
warning('error', 'MATLAB:nearlySingularMatrix');
warning('error', 'MATLAB:illConditionedMatrix');
warning('off', 'MATLAB:rankDeficientMatrix');

% determine the order and degree of associated spherical harmonics
[m, n] = sph_harm_ind_list(shOrder);
nParams = length(n);

%% Load variables for peak reconstruction

dataDir = fullfile(fileparts(mfilename('fullpath')), 'reconstruction_basis');

% load peak_indices reconstruction basis
reconBasis = dlmread(fullfile(dataDir, 'reconstruction_basis.txt'));
nPointsBasis = size(reconBasis, 1)/2;
reconBasis = reconBasis(1:nPointsBasis, :);
faces = dlmread(fullfile(dataDir, 'faces.txt'));
faces = faces - (faces > nPointsBasis)*nPointsBasis;

[~, reconBasis_theta, reconBasis_phi] = ...
    cart2sphere(reconBasis(:, 1), reconBasis(:, 2), reconBasis(:, 3));

reconBasis_sh = real_sph_harm(m, n, ...
    reconBasis_theta, reconBasis_phi)';


%% Load reference sphere for CSD algorithm
% load reference sphere to which we sample.

refSphere = dlmread(fullfile(fileparts(mfilename('fullpath')), ...
    'csd_data', 'sphere_small.txt'));

[~, ref_theta, ref_phi] = ...
    cart2sphere(refSphere(:, 1), refSphere(:, 2), refSphere(:, 3));

ref_sh = real_sph_harm(m, n, ref_theta , ref_phi)';


%% Prepare data

nWeightedScans = nnz(gtab.bvals > 0);
nVoxels = size(signal_intensities, 1);
nScans = size(gtab.bvals, 1);

if nonlinearitiesFlag
    % Correct expanded bvecs and bvalues (one set for each voxel).
    % TODO: vectorize for-loop for speed improvements.
    bvals = zeros(nScans, nVoxels);
    bvecs = zeros(nScans, 3, nVoxels);
    for iV = 1:nVoxels
        I = eye(3);
        v = gtab.bvecs*(I+reshape(nonlinearities(iV, :), [3 3]));
        vNorm = sqrt(sum(v.^2, 2));
        
        bvecs(:, :, iV) = v ./ vNorm;
        bvecs(vNorm == 0, :, iV) = 0;
        bvals(:, iV) = vNorm.^2.*gtab.bvals;
    end
    gtab.bvecs = reshape(permute(bvecs, [1 3 2]), [], 3);
    gtab.bvals = reshape(bvals, [], 1);
end

% convert gradient vectors to b-vectors
bvals = gtab.bvals;
bvecs = gtab.bvecs;
weightedScans = bvals > 0;

% convert cartesian bvec directions to spherical coordinates
[~, S_theta, S_phi] = cart2sphere(bvecs(weightedScans, 1), ...
    bvecs(weightedScans, 2), ...
    bvecs(weightedScans, 3));

% determine spherical harmonics of bvecs.
Q = real_sph_harm(m, n, S_theta , S_phi)';

% model
response_DW = estimate_response(bvecs, bvals, ...
    response.evals, response.S0);
response_sh = linsolve(Q, response_DW(weightedScans));

% Calculate the rotational harmonic decomposition up to
% harmonic order `m`, degree `n` for an axially and antipodally
% symmetric function. Note that all ``m != 0`` coefficients
% will be ignored as axial symmetry is assumed. Hence, there
% will be ``(sh_order/2 + 1)`` non-zero coefficients.
mask = m == 0;
dirac_sh = real_sph_harm(m(mask), n(mask), 0, 0);
response_rh = response_sh(mask) ./ dirac_sh;

% put r_rh in diagonal matrix.
R = diag(response_rh(floor(n ./ 2)+1));

% describe normalization lambda as in Dipy
% https://github.com/dipy/dipy/pull/340/files/6e8ef5fe90b01cfb048d36ef02f277b4029f165d#
lambda = (lambda * size(R, 1)  * response_rh(1) / ...
    (sqrt(size(ref_sh, 1)) * sqrt(362)));

B_ref_norm = lambda * ref_sh;
X = Q * R;

if nonlinearitiesFlag
    % reshape X to 3D tensor (scans X basis * voxels).
    X = reshape(X, nWeightedScans, nVoxels, []);
    Xall = permute(X, [1 3 2]);
    weightedScans = weightedScans(1:size(signal_intensities, 2));
end


%% computation

B_ref_norm = single(B_ref_norm);
peak_indices = zeros(size(signal_intensities, 1), outputPeaks, 'uint16');

if test.status
    test.shm_coeff = zeros(size(signal_intensities, 1), nParams);
end

exit_errors_all = categorical(repmat({''}, size(signal_intensities, 1), 1));
for i = 1:size(signal_intensities, 1)
    
    % signal intensities for the current voxel
    s = signal_intensities(i, :)';
    
    if nonlinearitiesFlag
       X = Xall(:, :, i);
    end
    
    % remove unsuccessful measurements
    if all(s(~weightedScans) <= 0)
        continue;
    end
    
    [shm_coeff, exit_error] = csdeconv(s(weightedScans), X, B_ref_norm, ...
        tau);
    exit_errors_all(i) = exit_error;
    
    if test.status
        test.shm_coeff(i, :) = shm_coeff;
    end
    
    odf = reconBasis_sh * shm_coeff;
    
    % find local maxima (adapted from
    % dsi-studio.labsolver.org/Manual/Reconstruction)
    peak_values = odf;
    peak_values(faces(odf(faces(:, 2)) >= odf(faces(:, 1)) | ...
        odf(faces(:, 3)) >= odf(faces(:, 1)), 1)) = 0;
    peak_values(faces(odf(faces(:, 1)) >= odf(faces(:, 2)) | ...
        odf(faces(:, 3)) >= odf(faces(:, 2)), 2)) = 0;
    peak_values(faces(odf(faces(:, 1)) >= odf(faces(:, 3)) | ...
        odf(faces(:, 2)) >= odf(faces(:, 3)), 3)) = 0;
    detected_peaks = peak_values./max(peak_values) > minPeakRatio;
    
    % skip isotropic voxels
    numberOfPeaks = sum(detected_peaks);
    if numberOfPeaks > maxPeaks
        continue;
    end
    
    % convert peak values into peak indices
    peak_values = peak_values .* detected_peaks;
    
    if any(peak_values)
        [~, I] = sort(peak_values, 'descend');
        m = min(outputPeaks, nnz(peak_values));
        peak_indices(i, :) = [I(1:m)' zeros(1, outputPeaks - m)];
    end
    
end   

% find peaks
reconBasis = [0 0 0; reconBasis];
diffusion_peaks = zeros(size(peak_indices, 1), 3, size(peak_indices, 2));
for iP = 1:size(peak_indices, 2)
    diffusion_peaks(:, :, iP) = reconBasis(peak_indices(:, iP)+1, :);
end

if any(~ismissing(exit_errors_all))
   fprintf('Warning:\n');
   exit_error_messages = categories(exit_errors_all);
   exit_errors_count = countcats(exit_errors_all);
   for i=1:length(exit_error_messages)
       fprintf('In %i voxels: %s\n', exit_errors_count(i), ...
           exit_error_messages{i});
   end
end

warning(warningsOld);


