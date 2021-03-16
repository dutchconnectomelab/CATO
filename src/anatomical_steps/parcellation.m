function parcellation(configParams)
% PARCELLATION   Execute parcellation scripts.
%
%   parcellation(CONFIGPARAMS) executes the parcellation script associated
%   with each template as specified in CONFIGPARAMS.

%% Initialization

status.parcellation = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---parcellation started----\n');

freesurferDir = configParams.general.freesurferDir;
templates = configParams.general.templates;

configParamsList = struct2list(configParams);
referenceFile = configParamsList{contains(configParamsList(:, 1), ...
    'ReferenceFile', 'IgnoreCase', true), 2};
registrationMatrixFile = configParamsList{contains(configParamsList(:, 1), ...
    'registrationMatrixFile', 'IgnoreCase', true), 2};

forceFreesurferOverwrite = jsonencode(configParams.parcellation.forceFreesurferOverwrite);

matchROIsFlag = configParams.parcellation.matchROIs;


%% Execute parcellation scripts for each template
for iTemplate = 1:length(templates)
    
    thisTemplate = templates{iTemplate};
    fprintf('template: %s\n', thisTemplate);
    
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
    
    exitCode = system(cmd);
    
    if exitCode ~=0
        error('CATO:parcellation:errorInParcellationScript', ...
            'Error in executing parcellation script:\n %s', ...
            cmd);
    end
    
    % Script creates:
    %     collect_region_properties.statsLhFile
    %     collect_region_properties.statsRhFile
    %     Surface annotation file freesurferDir/label/?h.TEMPLATE.annot
    %     Volume annotation file freesurferDir/mri/TEMPLATE+aseg.mgz
    %     ParcellationFile
    
    if matchROIsFlag
       
        fprintf('Match ROIs in parcellation file with color lookup table.\n');
        
        % Load brain parcellation (volume mapped to dwiReferenceFile)
        parcellation = load_nifti(parcellationFile);
        
        colorLookupTableFile = strrep(configParams.parcellation.lutFile, ...
            'TEMPLATE', thisTemplate);
        
        LUT = readColorLookupTable(colorLookupTableFile);
        
        % LEFT
        % Read original annotation file for orginal LUT
        annotFile = fullfile(freesurferDir, ...
            strrep('label/lh.TEMPLATE.annot', 'TEMPLATE', thisTemplate));
        [~, ~, oldLUT] = read_annotation(annotFile);
        
        % Match and reorder 
        [~, I] = ismember(strcat('ctx-lh-', oldLUT.struct_names), LUT.regionDescriptions);
        
        indxLeft = (parcellation.vol >=1000) & (parcellation.vol <2000);
        valuesLeft = parcellation.vol(indxLeft);
        parcellation.vol(indxLeft) = LUT.ROIs(I(valuesLeft - 1000 + 1));
        
        % RIGHT
        % Read original annotation file for orginal LUT
        annotFile = fullfile(freesurferDir, ...
            strrep('label/rh.TEMPLATE.annot', 'TEMPLATE', thisTemplate));
        [~, ~, oldLUT] = read_annotation(annotFile);
        
        % Match and reorder 
        [~, I] = ismember(strcat('ctx-rh-', oldLUT.struct_names), LUT.regionDescriptions);
        
        indxRight = (parcellation.vol >= 2000) & (parcellation.vol < 3000);
        valuesRight = parcellation.vol(indxRight);
        parcellation.vol(indxRight) = LUT.ROIs(I(valuesRight - 2000 + 1));   
        
        save_nifti(parcellation, parcellationFile);
        
    end
        
end

%% Clean up

status.parcellation = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---parcellation finished----\n');
