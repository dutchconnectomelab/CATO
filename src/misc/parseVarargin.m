function [configFile, runType, configParamsCl] = parseVarargin(varargin)
% PARSEVARARGIN    Parse input variables of structural and functional pipelines
%
%   INPUT VARIABLES
%   varargin:
%   Input variables of both pipelines.
%
%   OUTPUT VARIABLES
%   configFile:
%   Configuration file provided on command line (empty if not provided).
%
%   runType:
%   Run type provided on command line ('none' if not provided).
%
%   configParamsCl:
%   Additional name-value parameter pairs (in list-format).

assert(mod(numel(varargin), 2) == 0, ...
    'CATO:WrongNumberArgs', ...
    'Optional input arguments must come in pairs.');

% parse optional name-value pair arguments
configFile = [];
runType = 'none';
configParamsCl = cell(0,2); % NPARAMS x [variable NAME, variable VALUE]

while ~isempty(varargin)
    
    varName = varargin{1};
    varValue = varargin{2};
    
    switch lower(varName)
        
        case 'configurationfile'
            
            assert(ischar(varValue), ...
                'CATO:configurationFileNotText', ...
                ['configurationFile must be a row vector ', ...
                'of characters or string scalar.']);
            assert(isfile(varValue), ...
                'CATO:configurationFileNotFile', ...
                'configurationFile (%s) is not a file.', varValue);
            configFile = varValue;
            
        case 'runtype'
            
            varValue = lower(varValue);
            assert(ismember(runType, {'none', 'continue', 'overwrite'}));
            runType = varValue;
            
        case 'reconstructionsteps'
            
            % Command line uses ',' to split reconstruction steps.
            if ischar(varValue)
                varValue = strsplit(varValue, ',');
            end
            
            configParamsCl(end+1, :) = [{'general.reconstructionSteps'} {varValue}];
            
        otherwise
            
            % Command line uses -- to indicate variable names
            varName = erase(varName, '--');
            
            % Update varValue if it represents a numeric, logical or cell variable
            if ischar(varValue)
                try
                    varValueDecoded = jsondecode(varValue);
                    if isnumeric(varValueDecoded) || iscell(varValueDecoded) || islogical(varValueDecoded)
                        varValue = varValueDecoded;
                    end
                end
            end
            
            configParamsCl(end+1, :) = [varName {varValue}];
            
    end
    
    % remove the two entries we have just dealt with
    varargin(1:2) = [];
    
end
