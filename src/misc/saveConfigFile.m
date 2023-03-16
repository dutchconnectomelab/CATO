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

[fid, errorMessage] = fopen(configFile, 'w+');

if fid == -1
    error('CATO:saveConfigFile:cannotWriteToFile', 'Cannot write to configuration file (''%s''):\n%s.', configFile, errorMessage);
end

fprintf(fid, '%s', configParams);

fclose(fid);

end
