function ES = estimate_response(bvecs, bvals, evals, S0)
% ES = ESTIMATE_RESPONSE(bvecs, bvals, evals, S0) compute the expected 
% response ES (nweighted x 1) for each direction for a fiber in the 
% z-axis direction.
%
% INPUT:
% bvecs
%
% bvals
%
% evals
% Eigenvalues of the diffusion tensor. 
%
% S0
% Strength of signal in the presence of no diffusion gradient 
% (also called the b=0 value).
%
% OUTPUT
% ES
% Simulated signal: ES(q, tau) = S_0 e^(-b g^T R D R.T g).

% evecs describes the direction of the fiber (up along z-axis)
fiber_evecs = [0 0 1; 0 1 0; 1 0 0];

% now estimate the response of such fiber to each weighting.
n_gradients = size(bvecs, 1);
evals = evals(:);
RD = fiber_evecs * diag(evals) * fiber_evecs';

ES = nan(n_gradients, 1);

for i = 1:n_gradients
    g = bvecs(i, :)';
    ES(i) = S0 .* exp(-bvals(i) .* g' * RD * g);
end
