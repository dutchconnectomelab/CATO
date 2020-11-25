function peakList = tensor2peak(diffusionTensors)
% TENSOR2PEAK   Find peaks associated with tensors.
%
%   PEAKLIST = tensor2peak(DIFFUSIONTENSORS) calculates for all diffusion
%   tensors in DIFFUSIONTENSORS (N x 7) the associated peaks provided in
%   PEAKLIST (N x 3).

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


peakList = zeros(nTensors, 3);
for iTensor = 1:nTensors
    thisD = D(:, :, iTensor);

    if all(isnan(diag(thisD)))
        continue
    end
    
    [V, E] = eig(thisD);
    E = diag(E);
    [~, I] = sort(E, 'descend');
        
    peakList(iTensor, :) = V(:, I(1));
    
end

assert(all(isreal(peakList)), ['Diffusion peaks (peakList)', ...
    'contains non-real entries. Something is going wrong.']);

end
