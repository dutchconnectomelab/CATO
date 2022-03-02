function gtab = load_gtab(bvalsFile, bvecsFile, bValueZeroThreshold, bValueScalingTol)
warning off backtrace % TODO: better warning handling

% Load b-values file.
gtab.bvals = dlmread(bvalsFile);
validateattributes(gtab.bvals, ...
    {'numeric'}, {'nonempty', 'vector', 'nonnegative'}, ...
    mfilename, '(possibly processed) bvals');
gtab.bvals = gtab.bvals(:);

% Load b-vectors file.
gtab.bvecs = dlmread(bvecsFile);
validateattributes(gtab.bvecs, ...
    {'numeric'}, {'nonempty', '2d'}, ...
    mfilename, '(possibly processed) bvecs');

assert(any(size(gtab.bvecs) == 3), ...
    'B-vectors (%ix%i) must be a Nx3 or 3xN matrix.', ...
    size(gtab.bvecs, 1), size(gtab.bvecs, 2));

if size(gtab.bvecs, 1) == 3
    gtab.bvecs = gtab.bvecs';
end

if (min(gtab.bvals) > bValueZeroThreshold)
    warning(['No b0-scans identified. All b-values are higher than ', ...
        'the bValueZeroThreshold (%i) threshold. ', ...
        'Pipeline might not work.'], bValueZeroThreshold);
elseif any((gtab.bvals <= bValueZeroThreshold) & (gtab.bvals > 0))
    warning(['One or more b-values were below the bValueZeroThreshold', ...
        ' (%i) but nonzero. These are set to b-value = 0'], ...
        bValueZeroThreshold);
    
end

gtab.bvals(gtab.bvals <= bValueZeroThreshold) = 0;

% Check that bvectors are unit vectors (or zero).
inxWeighted = gtab.bvals > 0;
gtab.bvecs(~inxWeighted, :) = 0; % b0-scans no vector.

bvecsNorm = sqrt(sum(gtab.bvecs.^2,2));
if abs(1- bvecsNorm(inxWeighted)) >= bValueScalingTol
    warning(['b-vectors should be unit vectors. One or more vectors ', ...
        'deviates more from 1 than tolerance parameter bValueScalingTol']); 
end

% Check b-values and b-vectors match.
assert(isequal(length(gtab.bvals), size(gtab.bvecs, 1)), ...
    'CATO:load_gtab:bValsAndbVecsDoNotMatch', ...
    'Number of processed b-values (%i) and b-vectors (%i) do not match.', ...
    length(gtab.bvals), size(gtab.bvecs, 1));

end