function [fodf_sh, exitError] = csdeconv(dwi_data, A, B_reg, tau)
% [fodf_sh, i_num] = CSDECONV(dwi_data, A, B_reg, tau) perform the
% constrained-regularized spherical deconvolution.
%
% deconvolves the axially symmetric single fiber response function in
% rotational harmonic coefficients from the diffusion weighted signal.
%
% INPUT
% dwi_data:
% (NWEIGTHING x 1) diffusion weighted signals of a voxel to be deconvolved.
%
% A:
% (NWEIGHTING x NPARAMS) prediction matrix which estimates diffusion 
% weighted signals from FOD coefficients.
%
% B_reg:
% (NPOINTS x NPARAMS) SH basis matrix which maps FOD coefficients to FOD 
% values on the surface of the sphere. 
% B_reg should be scaled to account for lambda.
%
% tau:
% Threshold controlling the amplitude below which the corresponding fODF
% is assumed to be zero.  Ideally, tau should be set to zero. However, to
% improve the stability of the algorithm, tau is set to tau*100 % of the
% max fODF amplitude (here, 10% by default). This is similar to peak
% detection where peaks below 0.1 amplitude are usually considered noise
% peaks. Because SDT is based on a q-ball ODF deconvolution, and not
% signal deconvolution, using the max instead of mean (as in CSD), is
% more stable.
%
% OUTPUT
% fodf_sh:
% (NPARAMS x 1) spherical harmonics coefficients of the constrained-
% regularized fiber ODF.

%     https://github.com/dipy/dipy/blob/e34f71ced4956d6d4a1d364cfb57c42657e93e3e/dipy/reconst/csdeconv.py#L467


convergence = 50; % max number iterations
exitError = '';

% We do not use the precomputation trick, because this resulted in numerical
% errors when using a high number of harmonics (in sr-csd)
fodf_sh = A \ dwi_data;

% For the first iteration we use a smooth FOD that only uses SH orders up
% to 4 (the first 15 coefficients).
fodf = B_reg(:, 1:15) * fodf_sh(1:15);

% The mean of an fodf can be computed by taking $Y_{0,0} * coeff_{0,0}$
% tau was set to 10% of the mean FOD amplitude.
threshold = B_reg(1, 1) * fodf_sh(1) * tau;
threshold = max(threshold, 0);
where_fodf_small = find(fodf < threshold);

% If the low-order fodf does not have any values less than threshold, the
% full-order fodf is used.
if isempty(where_fodf_small)
    fodf = B_reg * fodf_sh;
    where_fodf_small = find(fodf < threshold);
    % If the fodf still has no values less than threshold, return the fodf.
    if isempty(where_fodf_small)
        i_num = 0;
        return
    end
end

for i_num = 1:convergence
    % This is the super-resolved trick.
    
    lambdaL = B_reg(where_fodf_small, :);
   
    M = [A; lambdaL];
    ba = [dwi_data; zeros(size(lambdaL, 1), 1)];

    
    if size(M,1) <= size(M, 2)
        break
    end
    
    fodf_sh = M \ ba;
    
    % Sample the FOD using the regularization sphere and compute k.
    fodf = B_reg * fodf_sh;
    where_fodf_small_last = where_fodf_small;
    where_fodf_small = find(fodf < threshold);
    
    if isequal(where_fodf_small, where_fodf_small_last)
        break
    end
   
end

if i_num == convergence
    exitError = sprintf('Spherical deconvolution (csdeconv.m) failed to converge.');
end
