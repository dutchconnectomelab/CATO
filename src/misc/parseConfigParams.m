function configParams = parseConfigParams(configParams)
% PARSECONFIGPARAMS Replace special configuration variables
%
%   configParams = parseConfigParams(configParams) updates the special
%   configuration variables (listed below).
%
%   INPUT VARIABLES
%   configParams:
%   Configuration parameters in list- or struct-format.
%
%   OUTPUT VARIABLES
%   configParams:
%   Configuration parameters in struct-format in which special variables
%   are replaced.
%
%   NOTE
%   Special variables that are replaced:
%   1. TOOLBOXDIR - replaced with toolbox directory
%   2. TEMPDIR - replaced with temporary directory location.
%   3. Variables in capitals (except FAMeasure) - replaced with variable value.
%   4. Reconstruction steps - transformed to cell-format.
%
%   Special variables must be seperated by a non-text character. For
%   example MAXFD_MAXDVARS.

% parse special variables
if isstruct(configParams)
    configParamsList = struct2list(configParams);
else
    configParamsList = configParams;
end

% replace TOOLBOXDIR
indx = strcmp(configParamsList(:, 1), 'general.TOOLBOXDIR');
if ~any(indx)
    TOOLBOXDIR = fileparts(fileparts(mfilename('fullpath')));
else
    TOOLBOXDIR = configParamsList{indx,2};
end
configParamsList(:, 2) = replacevariable(configParamsList(:, 2), ...
    'TOOLBOXDIR', TOOLBOXDIR);

% replace TEMPDIR (if defined).
indx = strcmp(configParamsList(:, 1), 'general.TEMPDIR');
if any(indx)
    if isempty(configParamsList{indx, 2})
        TEMPDIR = '/';
        while isfolder(TEMPDIR)
            TEMPDIR = fullfile(tempdir, ['CATO_' dec2hex(randi(2^20))]);
        end
        configParamsList{indx, 2} = TEMPDIR;
    end
end

% replace variables in capitals except FAMeasure.DTI FAMeasure.GQI etc.
for iP = 1:size(configParamsList, 1)
    specialVar = configParamsList{iP,1};
    
    if contains(specialVar, 'FAMeasure')
        continue
    end
    
    specialVar = strsplit(specialVar, '.');
    specialVar = specialVar{end};
    specialVar = upper(specialVar);
    configParamsList(:, 2) = replacevariable(configParamsList(:,2), ...
        specialVar, configParamsList{iP,2});
end

% replace make methods in small letters
indx = strcmp(configParamsList(:, 1), 'general.reconstructionMethods');
if any(indx)
    configParamsList{indx, 2} = lower(configParamsList{indx, 2});
end

% transform char into cell
indx = find(ismember(configParamsList(:, 1), ...
    {'general.reconstructionMethods', ...
    'general.templates', ...
    'reconstruction_diffusion.exportNifti.measures'}));
if any(indx)
    for i = 1:length(indx)
        if ~iscell(configParamsList{indx(i),2})
            configParamsList(indx(i), 2) = {{configParamsList{indx(i), 2}}};
        end
    end
end

% Reconstruction steps can be provided as char, struct or cell. 
configParams = list2struct(configParamsList);

if isfield(configParams, 'general')
    if isfield(configParams.general, 'reconstructionSteps')
    switch class(configParams.general.reconstructionSteps)
        case 'char'
            configParams.general.reconstructionSteps = ...
                {configParams.general.reconstructionSteps};
        case 'struct'
            % Form: XX.STEP1 = true; XX.STEP2 = false; etc.
            tmp = struct2list(configParams.general.reconstructionSteps);
            reconstructionSteps = tmp(:,1);
            
            assert(all(cellfun(@islogical, tmp(:,2))), ...
                'CATO:reconstructionStepsWrongFormat', ...
                ['Variable reconstructionSteps must be either ', ...
                'a char (''step1''), a cell ({''step1'',''step2''}), ', ...
                'or a struct (reconstructionSteps.step1 = true)']);
            
            reconstructionStepsFlag = [tmp{:,2}];
            reconstructionSteps = reconstructionSteps(reconstructionStepsFlag);
            configParams.general.reconstructionSteps = reconstructionSteps';
        case 'cell'
            % output format; do nothing.
        otherwise
            error('CATO:structural_pipeline:parameterWrong', ...
                ['Variable general.reconstructionSteps must be ', ...
                'a char or cell array, or struct.']);
    end

    end
end

end

function C = replacevariable(C, old, new)

assert(iscell(C));

if ~ischar(old)
    return
end

if ~ischar(new)
    new = jsonencode(new);
end

indx = cellfun(@ischar, C);
% C(indx) = cellfun(@(x) strrep(x, old, new), ...
%     C(indx), 'UniformOutput', false);

% Match only perfect matching strings TEMPLATE vs TEMPLATESDIR
%  lookbehind (?<=) and look ahead (?=)
% [^] Any character not contained within the brackets
% Or let it be the beginning (^) or end ($).

% Replace metacharacters that modify REPLACE in regexprep
new = strrep(new, '\', '\\');
new = strrep(new, '$', '\$');

C(indx) = cellfun(@(x) regexprep(x, ['((?<=[^A-Z])|^)' old '((?=[^A-Z]|$))'], new), ...
    C(indx), 'UniformOutput', false);
end
