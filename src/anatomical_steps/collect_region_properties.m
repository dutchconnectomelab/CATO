function collect_region_properties(configParams)
% COLLECT_REGION_PROPERTIES     Create regionPropertiesFile from
%                               FreeSurfer output.

%   collect_region_properties(CONFIGPARAMS) collates the statistics from
%   Freesurfer (volume, surface area and thickness) and measures obtained
%   from the parcellationFile (center of mass x, y and z and number of
%   voxels per region) in the regionPropertiesFile.

%% Initialization
status.collect_region_properties = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---collect_region_properties started----\n');

%% Loop over templates

templates = configParams.general.templates;

for iTemplate = 1:length(templates)
   
    thisTemplate = templates{iTemplate};
    
    fprintf('template: %s\n', thisTemplate);
    
    % Prepare file names
    regionPropertiesFile = strrep(configParams.collect_region_properties.regionPropertiesFile, ...
        'TEMPLATE', thisTemplate);
    statsLhFile = strrep(configParams.collect_region_properties.statsLhFile, ...
        'TEMPLATE', thisTemplate);
    statsRhFile = strrep(configParams.collect_region_properties.statsRhFile, ...
        'TEMPLATE', thisTemplate);
    statsSubFile = strrep(configParams.collect_region_properties.statsSubFile, ...
        'TEMPLATE', thisTemplate);
    lutFile = strrep(configParams.collect_region_properties.lutFile, ...
        'TEMPLATE', thisTemplate);
    ROIsFile = strrep(configParams.general.ROIsFile, ...
        'TEMPLATE', thisTemplate);
    parcellationFile = strrep(configParams.parcellation.parcellationFile, ...
        'TEMPLATE', thisTemplate);
    
    statsFileReadOptions = {'CommentStyle', '#', ...
                            'filetype', 'text', ...
                            'ReadVariableNames', false};
    
    % Load left hemisphere ananatomical statistics
    Clh = readtable(statsLhFile, statsFileReadOptions{:});
    
    Clh.Properties.VariableNames = readColHeaders(statsLhFile);
    Clh.regionDescriptions = strcat('ctx-lh-', Clh.StructName);
    Clh.Properties.VariableNames{'GrayVol'} = 'volume';
    Clh.Properties.VariableNames{'SurfArea'} = 'surface_area';    
    Clh.Properties.VariableNames{'ThickAvg'} = 'thickness';    
    
    % Load left hemisphere ananatomical statistics
    Crh = readtable(statsRhFile, statsFileReadOptions{:});
    
    Crh.Properties.VariableNames = readColHeaders(statsRhFile);
    Crh.regionDescriptions = strcat('ctx-rh-', Crh.StructName);
    Crh.Properties.VariableNames{'GrayVol'} = 'volume';
    Crh.Properties.VariableNames{'SurfArea'} = 'surface_area';   
    Crh.Properties.VariableNames{'ThickAvg'} = 'thickness';    
    
    
    % Load subcortical anatomical statistics
    Csub = readtable(statsSubFile, statsFileReadOptions{:});
    
    varNames = readColHeaders(statsSubFile);
    Csub.Properties.VariableNames = [varNames; {'uknown'}];
    Csub.Properties.VariableNames{'StructName'} = 'regionDescriptions';
    Csub.Properties.VariableNames{'Volume_mm3'} = 'volume';
    Csub.thickness = nan(height(Csub), 1);
    Csub.surface_area = nan(height(Csub), 1);
    
    % Combine left, right and subcortical statistics
    propertiesOI = {'volume', 'surface_area', 'thickness', 'regionDescriptions'};
    C = [Clh(:, propertiesOI); Crh(:, propertiesOI); Csub(:, propertiesOI)];
    
    %% Merge statistics with ROIs lookup tbale
    LUT = readtable(lutFile,  'filetype', 'text', 'ReadVariableNames', false);
    LUT.Properties.VariableNames = {'ROIs', 'regionDescriptions', 'Color1', 'Color2', 'Color3', 'Other'};
    LUT.regionDescriptions = deblank(LUT.regionDescriptions);
    C = outerjoin(C, LUT, 'MergeKeys', true, 'Type', 'right');
    
    % Load ROIs file
    ROIs = dlmread(ROIsFile);
            
    %% Calculate center of mass statistics and number of voxels 
    parcellation = load_nifti(parcellationFile);

    % Get vox2ras-tkr matrix.
    % Use vox2ras matrix to determine the orientation of the voxels.
    % vox2ras-tkr is voxel*voxel size with center of all points as origin
    [~, J] = max(abs(parcellation.vox2ras(1:3, 1:3)), [], 2);
    voxelSize = parcellation.pixdim(2:4);
    dimRAS = voxelSize .* size(parcellation.vol)';

    T = zeros(4);
    T(1, J(1)) = sign(parcellation.vox2ras(1, J(1)))*voxelSize(1);
    T(2, J(2)) = sign(parcellation.vox2ras(2, J(2)))*voxelSize(2);
    T(3, J(3)) = sign(parcellation.vox2ras(3, J(3)))*voxelSize(3);
    T(1, 4) = -sign(parcellation.vox2ras(1, J(1)))*dimRAS(J(1))/2;
    T(2, 4) = -sign(parcellation.vox2ras(2, J(2)))*dimRAS(J(2))/2;
    T(3, 4) = -sign(parcellation.vox2ras(3, J(3)))*dimRAS(J(3))/2;
    T(4, 4) = 1;

    parcellation = parcellation.vol;
    
    coordinatesCenter = nan(length(ROIs), 3);
    nVoxels = nan(length(ROIs), 1);
    for ir = 1:length(ROIs)
        Vindices = find(parcellation == ROIs(ir));
        
        if ~isempty(Vindices)
            [Vc, Vr, Vs] = ind2sub(size(parcellation), Vindices);
            V = [Vc Vr Vs] - 1; % indexings starts at 0 in NIFTI
            coordinates = T*[V ones(size(V, 1), 1)]';
            coordinatesCenter(ir, :) = mean(coordinates(1:3, :), 2);
            nVoxels(ir) = numel(Vindices);
        end
        
    end
        
    XYZN = array2table([ROIs coordinatesCenter nVoxels], 'VariableNames', ...
        {'ROIs', 'coordinate_x', 'coordinate_y', 'coordinate_z', 'number_of_voxels'});

    % The x-, y-, z-coordinates are NOT centered at (0,0,0).
    
    C = outerjoin(C, XYZN, 'MergeKeys', true, 'Type', 'left');
   
    %% Combine and save
    
    % Ensure same properties and order as in previous versions.
    propertyDescriptions = {'coordinate_x', 'coordinate_y', ...
                            'coordinate_z', 'volume', ...
                            'surface area', 'thickness', ...
                            'number of voxels'}';
                        
    % save regionProperties
    [~, indxROIs] = ismember(ROIs, C.ROIs);
    [~, indxPropertyDescriptions] = ismember(...
        strrep(propertyDescriptions, ' ', '_'), C.Properties.VariableNames);
    
    if any(indxROIs == 0)
        error('ROI not included in template');
    end
    
    regionDescriptions = C.regionDescriptions(indxROIs);
    regionProperties = C{indxROIs, indxPropertyDescriptions}; %#ok
    
    regionPropertiesTable = C(indxROIs, indxPropertyDescriptions);
    regionPropertiesTable.Properties.RowNames = regionDescriptions;

    save(regionPropertiesFile, 'ROIs', 'propertyDescriptions', ...
                               'regionDescriptions', 'regionProperties', ...
                               'regionPropertiesTable');
                           
    
end

%% Clean up
status.collect_region_properties = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---collect_region_properties finished----\n');

end

function varNames = readColHeaders(statsFile)

    varNames = splitlines(fileread(statsFile));
    varNames = varNames(cellfun(@(x) contains(x, 'ColHeaders'), varNames));
    
    assert(numel(varNames) == 1, 'CATO:collect_region_properties', ...
        'Error reading stats file (%s).\n More than 1 header line.', statsFile);
    
    varNames = strtrim(erase(varNames, '# ColHeaders'));
    varNames = split(varNames{1});
    
end

