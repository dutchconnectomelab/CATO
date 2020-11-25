function response = auto_response(diffusionMeasures, weightDescriptions, ...
    S, bvals, CC_mask, FA_threshold)
% response = AUTO_RESPONSE(diffusion_measure_data, ...
%     signal_intensities, gradient_table, wm_mask) calcules the S0 and
%     eigenvalues of the response function.
%
% We base the response function on voxels with FA-values higher than 
% FA_threshold in the corpus collosum (as indicated in CC_mask).
%
% INPUT
% diffusion_measure_data
% (NX x NY x NZ x NMEASURES array) with  FA, RD, AD for every voxel.
%
% signal_intensities
% (NX x NY x NZ x NWEIGHTS array) with signal intensities.
%
% gradient_table
% (NWEIGHTS x 3) table with gradients.
%
% CC_mask
% (NX x NY x NZ array) with ones for corpus collosum voxels and zeros
% otherwise.
%
% OUTPUT
% response:
% Structure with 1) response.S0 describing the average intensity in b0
% scans and 2) response.evals a 3x1 vector with the first eigenvalue and
% average of second and third eigenvalue.



% load diffusion measures
FAindx = ismember(weightDescriptions, 'fractional anisotropy');
FA = diffusionMeasures(:, :, :, FAindx);

RDindx = ismember(weightDescriptions, 'radial diffusivity');
RD = diffusionMeasures(:, :, :, RDindx);

ADindx = ismember(weightDescriptions, 'axial diffusivity');
AD = diffusionMeasures(:, :, :, ADindx);

% apply FA threshold to signal intensities
indxFA = FA(:) > FA_threshold;
indxCC = CC_mask(:) > 0;

roi = indxFA & indxCC;

if nnz(indxFA) == 0
    error(['connectomizer:auto_response:' ...
        'NoVoxelsAboveFAThreshold'], ...
        ['No voxel with a FA higher than %2.2f were found.\n', ...
        'Try a larger roi or a lower threshold.'], FA_threshold);
end
if nnz(roi) == 0
    error(['connectomizer:auto_response:' ...
        'NoVoxelsAboveFAThreshold'], ...
        ['No voxel in the corpus callosum CCRegions mask found ', ...
        'with a FA higher than %2.2f were found.\n', ...
        'Try a larger roi or a lower threshold.'], FA_threshold);
end

% average signal intensity over nonweighted scans
nonweighted = bvals(:) < 0.1;

% obtain averaged values
S0 = mean(mean(S(roi, nonweighted)));

evals(1) = mean(AD(roi));
evals([2 3]) = mean(RD(roi));

% save in response variable
response.S0 = S0;
response.evals = evals;
