function saveConfigFile(configFile, configParams)
% SAVECONFIGFILE Write configuration parameters to file.
%
%   INPUT VARIABLES
%   configFile:
%   Configuration parameters are saved in this file.
%   
%   configParams:
%   Configuration parameters in struct-format.

configParams = jsonencode(configParams);

configParams = prettyjson(configParams);

fid = fopen(configFile, 'w+');

fprintf(fid, '%s', configParams);

fclose(fid);

end
