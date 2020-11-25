function configParams = readConfigFile(configFile)
% READCONFIGURATIONFILE Read configuration file.
%
%   INPUT VARIABLES
%   configFile:
%   Configuration file.
%
%   OUTPUT VARIABLES
%   configParams:
%   Configuration parameters in struct-format.

configParams = fileread(configFile);

try
configParams = jsondecode(configParams);
catch ME
    error('CATO:readConfigFile:JSONerror', ...
        'Cannot parse json syntax of file %s.\n%s', configFile, ME.message);
end

end

