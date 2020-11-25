function parcellation(configParams)
% PARCELLATION   Execute parcellation scripts.
%
%   parcellation(CONFIGPARAMS) executes the parcellation script associated
%   with each template as specified in CONFIGPARAMS.

%% Initialization

status.parcellation = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---Parcellation started----\n');

freesurferDir = configParams.general.freesurferDir;
templates = configParams.general.templates;

configParamsList = struct2list(configParams);
referenceFile = configParamsList{contains(configParamsList(:, 1), ...
    'ReferenceFile', 'IgnoreCase', true), 2};
registrationMatrixFile = configParamsList{contains(configParamsList(:, 1), ...
    'registrationMatrixFile', 'IgnoreCase', true), 2};

forceFreesurferOverwrite = jsonencode(configParams.parcellation.forceFreesurferOverwrite);

%% Execute parcellation scripts for each template
for iTemplate = 1:length(templates)
    
    thisTemplate = templates{iTemplate};
    
    templateScript = strrep(configParams.parcellation.templateScript, ...
        'TEMPLATE', thisTemplate);
    
    parcellationFile = strrep(configParams.parcellation.parcellationFile, ...
        'TEMPLATE', thisTemplate);
    
    cmd = ['"', templateScript, '"', ...
        ' --freesurferDir="', freesurferDir, '"', ...
        ' --referenceFile="', referenceFile, '"', ...
        ' --registrationMatrixFile="', registrationMatrixFile, '"', ...
        ' --parcellationFile="', parcellationFile, '"', ...
        ' --forceFreesurferOverwrite="', forceFreesurferOverwrite, '"'];
    
    fprintf('Template: %s\n', thisTemplate);

    exitCode = system(cmd);
    
    if exitCode ~=0
        error('CATO:parcellation:errorInParcellationScript', ...
            'Error in executing parcellation script:\n %s', ...
            cmd);
    end
    
    % Script creates:
    %     collect_region_properties.statsLhFile
    %     collect_region_properties.statsRhFile
        
end

%% Clean up

status.parcellation = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---Parcellation finished----\n');
