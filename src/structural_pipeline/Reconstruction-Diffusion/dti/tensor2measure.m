function [diffusionMeasures, weightDescriptions] = tensor2measure(diffusionTensors)
% TENSOR2MEASURE    Compute diffusion measures associated with a tensor
%
%   tensor2measure(DIFFUSIONTENSORS) calculates fractional anisotropy,
%   axial diffusivity, radial diffusivity and mean diffusivity of diffusion
%   tensors in DIFFUSIONTENSORS
%
%   INPUT VARIABLES
%   diffusionTensors:
%   nVoxels x 7 matrix with the estimated diffusion tensor parameters.
%
%   OUTPUT VARIABLES
%   diffusionMeasures:
%   nVoxels x 4 matrix with calcualted diffusion measures.
%
%   weightDescriptions:
%   Cell with names of the calculated diffusion measures.

nTensors = size(diffusionTensors, 1);

D = zeros(3,3, nTensors);
D(1,1,:) = diffusionTensors(:, 1);
D(2,2,:) = diffusionTensors(:, 2);
D(3,3,:) = diffusionTensors(:, 3);
D(1,2,:) = diffusionTensors(:, 4);
D(2,1,:) = diffusionTensors(:, 4);
D(1,3,:) = diffusionTensors(:, 5);
D(3,1,:) = diffusionTensors(:, 5);
D(2,3,:) = diffusionTensors(:, 6);
D(3,2,:) = diffusionTensors(:, 6);

diffusionMeasures = zeros(nTensors, 4);
weightDescriptions = cell(4,1);
for iTensor = 1:nTensors
    thisD = D(:, :, iTensor);
    
    if all(isnan(diag(thisD)))
        continue
    end
    
    [~, E] = eig(thisD);
    E = diag(E);
    E = sort(E, 'descend');
    
    E(E<0) = 0;
    
    if all(E == 0)
        continue
    end
    
    % FA:
    weightDescriptions{1} = 'fractional anisotropy';    
    diffusionMeasures(iTensor, 1) = sqrt(1/2) * ...
        sqrt((E(1) - E(2))^2 + (E(2) - E(3))^2 + (E(1) - E(3))^2) / ...
        sqrt(E(1)^2 + E(2)^2+E(3)^2);
    
    % AD: largest eigenvalue
    weightDescriptions{2} = 'axial diffusivity';    
    diffusionMeasures(iTensor, 2) = E(1);
    
    % RD: average of two smallest eigenvalues
    weightDescriptions{3} = 'radial diffusivity';
    diffusionMeasures(iTensor, 3) = mean(E(2:3));
    
    % MD: average over all three eigenvalues
    weightDescriptions{4} = 'mean diffusivity';
    diffusionMeasures(iTensor, 4) = mean(E);
    
end


end
