function orientationStr = nifti_orientation(vox2ras)
% NIFTI_ORIENTATION computes the orientation string from the vox2ras
% matrix.
%
% Input:
%   vox2ras: a 4x4 matrix (typically the vox2ras parameter from the
%   NIFTI file header).
%
% Output:
%   orientation_str: a string indicating the orientation (e.g. 'RAS' or
%   'LPI').

% Notes:
% For more information on the vox2ras transformation see the FreeSurfer
% documentation:
% https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems

% Extract the direction cosine matrix (3x3)
R = vox2ras(1:3, 1:3);

% Compute the orientation codes
orientationCode = zeros(1, 3);
for i = 1:3
    [~, maxIndx] = max(abs(R(:,i)));
    if R(maxIndx, i) > 0
        orientationCode(i) = maxIndx;
    else
        orientationCode(i) = -maxIndx;
    end
end

% Convert orientation codes to orientation string
orientationStr = '???';
for i = 1:3
    switch orientationCode(i)
        case 1
            orientationStr(i) = 'R';
        case -1
            orientationStr(i) = 'L';
        case 2
            orientationStr(i) = 'A';
        case -2
            orientationStr(i) = 'P';
        case 3
            orientationStr(i) = 'S';
        case -3
            orientationStr(i) = 'I';
    end
end
end