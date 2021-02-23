function colorLookupTable = readColorLookupTable(colorLookupTableFile)
% READCONFIGURATIONFILE Read configuration file.
%
%   INPUT VARIABLES
%   colorLookupTableFile:
%   Color lookup table file.
%
%   OUTPUT VARIABLES
%   colorLookupTable:
%   Color lookup table with columns: ROIs, regionDescriptions, R, G, B, A


% Load look up table.
fid = fopen(colorLookupTableFile);

try
    colorLookupTable = textscan(fid, '%d %s %d %d %d %d', inf, ...
        'CommentStyle', '#', ...
        'ReturnOnError', false);
    assert(numel(unique(cellfun(@length, colorLookupTable))) == 1, ...
        'Length of items in color lookup table varies.');
catch ME
    error('CATO:readColorLookupTable:parsingError', ...
        'Cannot parse color lookup table %s.\n%s', colorLookupTableFile, ME.message);
end

fclose(fid);

colorLookupTable = table(colorLookupTable{:}, ...
    'variableNames', {'ROIs', 'regionDescriptions', 'R', 'G', 'B', 'A'});

end