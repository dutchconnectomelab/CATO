function writeFibers(fibers, fiberFile, voxelSize, header)
% WRITEFIBERS   Write fibers to TRK-file
%
%   INPUT VARIABLES
%   fibers:
%   Variable with fibers in array format (nVoxels x 3 x nFibers) or cell
%   format (nFibers x 1).
%
%   fiberFile:
%   TRK file in which the reconstructed fiber cloud is saved. See:
%   http://trackvis.org/docs/?subsect=fileformat
%
%   voxelSize:
%   Voxel size of dwiProcessedFile
%
%   header (optional):
%   If the header for the fiberFile is provided, then a new fiberFile is
%   created starting with the header. If the header is not provided, then
%   the fiberFile is appended (expecting that the fiberFile with header is
%   created earlier).

%% Initialization

% If header is provided then create new fiberFile and write header to file.
if nargin == 4
    header.n_properties = 0;
    header.n_scalars = 0;
    
    [fid, ~] = fopen(fiberFile, 'Wb');
    writeTrkHeader(fid, header);
else
    [fid, ~] = fopen(fiberFile, 'Ab');
end

if isempty(fibers)
    fclose(fid);
    return
end

% Find the number of fibers.
if iscell(fibers)
    nFibers = length(fibers);
else
    nFibers = size(fibers, 3);
end

%% Write fibers to file.
for iFiber = 1:nFibers
    
    if iscell(fibers)
        
        thisFiber = fibers{iFiber};
        points = size(thisFiber, 2);
        
    else
        
        % Remove the zeros padding from the fiber.
        % SPEED: this is done when writing to optimize.
        %     thisFiber = thisFiber(:,thisFiber(1,:) ~=0);
        
        indx = fibers(1, :, iFiber) ~=0;
        points = nnz(indx);
        
        % Voxel space to coordinates
        % Trackvis uses the voxmm format. Points are saved by multiplying the voxel
        % coordinate with the voxel sizes. See for more information:
        % https://github.com/nipy/nibabel/blob/
        % 65d5fc61545f55a50a45a07fbbaeb99c2dbe6bbb/nibabel/trackvis.py#L343
        thisFiber = (fibers(:, indx, iFiber) -1) .* voxelSize;
        
    end
    
    % write reconstructed fiber to file
    fwrite(fid, points, 'int');
    fwrite(fid, thisFiber, 'float');
    
end


%% Close fiber file.
fclose(fid);

