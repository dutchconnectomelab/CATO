function [fibers, header] = readTrk(fiberFile)
% READTRK Read fiber file (in TRK-format).
%
%   INPUT VARIABLES
%   fiberFile:
%   Fiber cloud file (in TRK-format).
%
%   OUTPUT VARIABLES
%   fibers:
%   Fibers in cell-format.
%
%   header:
%   Header of fiber file.

[fid, message] = fopen(fiberFile, 'r');
header = readTrkHeader(fid);

iFiber = 0;

fibers = cell(1e7,1);

while true
    
    % read fiber
    iFiber = iFiber + 1;
    points = fread(fid, 1, 'int');
    thisFiber = fread(fid, [(3 + header.n_scalars) points], 'float');
    
    if header.n_properties > 0
        fread(fid, header.n_properties, 'float');
    end    
    
    if isempty(points)
        break;
    end

    fibers{iFiber} = thisFiber;
    
end

fclose(fid);

fibers = fibers(1:iFiber-1);