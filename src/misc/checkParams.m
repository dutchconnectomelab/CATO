function configError = checkParams(configParams, reconSteps, parameterPropertiesFile)
% CHECKPARAMS   Validate that all parameters confirm to the expected format
%
%   INPUT VARIABLES
%   configParams:
%   Structure with configuration parameters.
%
%   reconSteps:
%   Steps of which the configuration parameters are checked.
%
%   parameterPropertiesFile:
%   Properties file (.xlsx) that describes parameter names, attributes and
%   types.
%
%   OUTPUT VARIABLES
%   configError:
%   Cell with errors encountered during parameter validation.
%
%   NOTES
%   checkParams has a few special cases:
%   1. reconstructionSteps parameter
%      The pipeline checks that the steps correspond to a list.
%   2. FAMeasure
%      FAMeasure is coded in the configuration file as a struct, therefore
%      a 'general.FAMeasure' variable is derived.
%
%   and uses 2 additional attributes that extend 'validateattributes':
%   1. isFolder
%   2. isFile

assert(isstruct(configParams))
assert(ischar(reconSteps) || iscellstr(reconSteps) || isstring(reconSteps) || iscell(reconSteps))

if ~iscell(reconSteps)
    reconSteps = {reconSteps};
end

paramProps = readtable(parameterPropertiesFile);
configParamsList = struct2list(configParams);

configError = cell(size(configParamsList,1),4);
configError(:,1) = configParamsList(:,1);
configError(:,4) = cellfun(@jsonencode, configParamsList(:,2), 'uniformOutput', false);

% Functional pipeline remove double parameters
configParamsList(:, 1) = cellfun(@(x) regexprep(x, '_[0-9].', '.'), ...
    configParamsList(:, 1), 'UniformOutput', false);

for i = 1:size(configParamsList,1)
    
    % Special case: reconstructionSteps parameter
    if contains(configParamsList(i,1), 'reconstructionSteps')
        reconstructionSteps = configParamsList{i,2};
        if ~all(ismember(reconstructionSteps, ...
                {'structural_preprocessing', 'parcellation', ...
                'collect_region_properties', 'reconstruction_diffusion', ...
                'reconstruction_fibers', 'reconstruction_fiber_properties', ...
                'reconstruction_network', ...
                'functional_preprocessing', 'parcellation', ...
                'collect_region_properties', 'compute_motion_metrics', ...
                'reconstruction_functional_network'}))
            warning(['general.reconstructionSteps parameter (%s) ', ...
                'contains unrecognized reconstruction steps.'], ...
                jsonencode(reconstructionSteps));
        end
        for iR = 1:length(reconstructionSteps)
           assert(exist(reconstructionSteps{iR}, 'file') == 2, ...
               ['reconstructionStep (%s) is not a function in the ', ...
               'current folder or on the search path.'], reconstructionSteps{iR});    
        end
        continue
    end
    
    % Special case:  FAMeasure
    if contains(configParamsList(i,1), 'FAMeasure')
        % rename specific methods (general.FAMeasure.DTI)
        configParamsList{i,1} = 'general.FAMeasure';
    end
    
    indx = strcmp(configParamsList(i,1), paramProps.name);
    
    if all(~indx)
        configError{i, 3} = ['Variable properties not defined in ', ...
            'parameter properties file.'];
        continue
    end
    
    thisParamVal = configParamsList{i,2};
    
    % ValidateAttributes
    thisAttributes = paramProps.attributes{indx};
    thisAttributes = erase(thisAttributes, {'isfile', 'isfolder'});
    thisAttributes = strsplit(strtrim(thisAttributes));    
    numInAttributes = cellfun(@str2double, thisAttributes, 'uniformOutput', false);
    numInAttributesIndx = ~cellfun(@isnan, numInAttributes);
    thisAttributes(numInAttributesIndx) = numInAttributes(numInAttributesIndx);
    
    thisAttributes = thisAttributes(~cellfun(@isempty, thisAttributes));
    
    thisType = paramProps.type{indx};
    thisType = strsplit(strtrim(thisType));
    
    inputParamFlag = any(contains(paramProps.inputFile{indx}, reconSteps));
    outputParamFlag = any(contains(paramProps.outputFile{indx}, reconSteps));
    parameterParamFlag = any(contains(paramProps.parameter{indx}, reconSteps));
    
    if inputParamFlag || outputParamFlag || parameterParamFlag
        try
            validateattributes(thisParamVal, ...
                thisType, ...
                thisAttributes);
        catch ME
            configError{i, 2} = strjoin(splitlines(ME.message));
        end
    end
    
    % Validate additional attributes (isfile)
    % File is not generated by pipeline
    if inputParamFlag && ~outputParamFlag
        
        paramValExp = expandMethodFiles(thisParamVal, ...
            configParams.general.reconstructionMethods);
        paramValExp = expandTemplateFiles(paramValExp, ...
            configParams.general.templates);
        
        paramValExp = paramValExp(:);
        
        for iParamValExp = 1:size(paramValExp, 1)
            
            thisParamValExp = paramValExp{iParamValExp};
            
            if contains(paramProps.attributes(indx), 'isfile')
                if ~isfile(thisParamValExp)
                    if isempty(configError{i,2})
                        configError{i,2} = 'File not found.';
                    else
                        configError{i,2} = [configError{i,2}, ' ', thisParamValExp];
                    end
                end
            end
            if contains(paramProps.attributes(indx), 'isfolder')
                if ~isfolder(thisParamValExp)
                    if isempty(configError{i,2})
                        configError{i,2} = 'Folder not found.';
                    else
                        configError{i,2} = [configError{i,2}, ' ', thisParamValExp];
                    end
                end
            end
        end
    end
end
end

function varNameExpanded = expandMethodFiles(varName, reconstructionMethods)

if ~ischar(varName)
    varNameExpanded = [];
    return
end

if contains(varName, 'METHOD')
    nMethods = length(reconstructionMethods);
    varNameExpanded = cell(nMethods, 1);
    for iMethod = 1:nMethods
        thisMethod = reconstructionMethods{iMethod};
        thisVarNameMethod = strrep(varName, 'METHOD', thisMethod);
        varNameExpanded{iMethod} = thisVarNameMethod;
    end
else
    varNameExpanded = {varName};
end
end

function varValExp = expandTemplateFiles(varValues, templates)

if ~iscellstr(varValues)
    varValExp = [];
    return
end

if contains(varValues, 'TEMPLATE')
    varValExp = cell(size(varValues,1), length(templates));
    for i=1:size(varValues, 1)
        thisVarValue = varValues{i};
        
        for iTemplate = 1:length(templates)
            thisTemplate = templates{iTemplate};
            thisVarValTemplate = strrep(thisVarValue, 'TEMPLATE', thisTemplate);
            varValExp{i, iTemplate} = thisVarValTemplate;
        end
        
    end
else
    varValExp = varValues;
    
end
end

