function reconstruction_network(configParams)
% RECONSTRUCTION_NETWORK   Reconstruct the structural connectivity
%                          matrices from the fiber properties.
%
%   reconstruction_network(CONFIGPARAMS) reconstruct the structural
%   connectivity matrices for each template and each method from the
%   fiberPropertiesFile according to parameters specified in CONFIGPARAMS.

%% Initialization
status.reconstruction_network = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_network started----\n');

reconstructionMethods = lower(configParams.general.reconstructionMethods);
FAMeasure = configParams.general.FAMeasure;
regionPropertiesFile = configParams.collect_region_properties.regionPropertiesFile;
fiberPropertiesFile = configParams.reconstruction_fiber_properties.fiberPropertiesFile;
maxAngleRad = deg2rad(configParams.reconstruction_network.maxAngleDeg);
minFA = configParams.reconstruction_network.minFA;
minLengthMM = configParams.reconstruction_network.minLengthMM;
connectivityMatrixFile = configParams.reconstruction_network.connectivityMatrixFile;

templates = configParams.general.templates;

%% Reconstruct networks for each template and each method.
for iTemplate = 1:length(templates)
    for iMethod = 1:length(reconstructionMethods)
        thisMethod = reconstructionMethods{iMethod};
        thisTemplate = templates{iTemplate};
        
        %% Initialize
        thisRegionPropertiesFile = strrep(regionPropertiesFile, ...
            'TEMPLATE', thisTemplate);
        thisFiberPropertiesFile =  strrep(strrep(fiberPropertiesFile, ...
            'METHOD', thisMethod), ...
            'TEMPLATE', thisTemplate);
        thisConnectivityMatrixFile = strrep(strrep(connectivityMatrixFile, ...
            'METHOD', thisMethod), ...
            'TEMPLATE', thisTemplate);
                
        % regionDescriptions
        data = load(thisRegionPropertiesFile);
        regionDescriptions = data.regionDescriptions;
        ROIs = data.ROIs;
        nRegions = length(regionDescriptions);
        
        indx = strcmp(data.propertyDescriptions, 'volume');
        regionVolumes = data.regionProperties(:, indx);
        regionVolumesAvg = repmat(regionVolumes, [1 nRegions]);
        regionVolumesAvg = 0.5*(regionVolumesAvg + regionVolumesAvg');           
        
        indx = strcmp(data.propertyDescriptions, 'surface area');        
        regionSurfaces = data.regionProperties(:, indx);
        regionSurfacesAvg = repmat(regionSurfaces, [1 nRegions]);
        regionSurfacesAvg = 0.5*(regionSurfacesAvg + regionSurfacesAvg');
        
        % Fiber properties
        data = load(thisFiberPropertiesFile);
        fiberProperties = data.fiberProperties;
        propertyDescriptions = data.propertyDescriptions;
        NPROPFIXED = 5;
        diffusionWeightDescriptions = propertyDescriptions(NPROPFIXED+1:end);
        nWeights = length(diffusionWeightDescriptions);
        
        
        %% Filter fibers
        filterIndx  = true(size(fiberProperties,1),1);
        
        % Maximum angle
        indx = strcmp(propertyDescriptions, 'maxAngle');
        filterIndx = filterIndx & (fiberProperties(:, indx) < maxAngleRad);

        % Minimum fiber length
        indx = strcmp(propertyDescriptions, 'length');
        filterIndx = filterIndx & (fiberProperties(:, indx) > minLengthMM);
        
        % Minimum FA measure
        % Dependent on the reconstruction method FA can be measured as FA or e.g.
        % generalized fractional anisotropy.
        FAMeasureMethods = fieldnames(FAMeasure);
        indxFAMeasure = strcmpi(FAMeasureMethods, thisMethod);
        if nnz(indxFAMeasure) == 0
            error(['Cannot find reconstruction method (%s) in ', ...
                'FAMeasure parameter (%s).'], ...
                thisMethod, lower(jsonencode(FAMeasure)));
        elseif nnz(indxFAMeasure) > 2
            error(['Multiple method (%s) instances found in ', ...
                'FAMeasure parameter (%s).'], ...
                thisMethod, lower(jsonencode(FAMeasure)));
        end
        
        indx = strcmpi(propertyDescriptions, ...
            FAMeasure.(FAMeasureMethods{indxFAMeasure}));
        filterIndx = filterIndx & (fiberProperties(:, indx) > minFA);        

        %% Construct matrices
        W = zeros(nRegions, nRegions, 1+nWeights+2);
        
        NOS = accumarray(fiberProperties(:, 2), double(filterIndx), ...
            [0.5*nRegions*(nRegions-1) 1]);
        % TODO: error if regionProperties file does not match
        % fiberProperties file.
        W(:, :, 1) = squareform(NOS);
        weightDescriptions{1} = 'number of streamlines';
        
        for iW = 1:nWeights
            thisW = accumarray(fiberProperties(:, 2), ... % indices
                filterIndx .* fiberProperties(:, NPROPFIXED+iW), ... % weightings
                [0.5*nRegions*(nRegions-1) 1]); % total size
            thisW = thisW ./ NOS;
            thisW(NOS == 0) = 0;
            
            W(:,:, 1+iW) = squareform(thisW);
            weightDescriptions(1+iW) = diffusionWeightDescriptions(iW);
        end
        
        % Add special measures
        % Streamline volume density
        W(:, :, 1+ nWeights + 1) = W(:, :, 1) ./ regionVolumesAvg;
        weightDescriptions{1 + nWeights + 1} = 'streamline volume density';
        
        % Streamline surface density
        W(:, :, 1+ nWeights + 2) = W(:, :, 1) ./ regionSurfacesAvg;
        weightDescriptions{1 + nWeights + 2} = 'streamline surface density';
        
        % Save (in correct format)
        connectivity = W;
        weightDescriptions = weightDescriptions(:); 
        
        save(thisConnectivityMatrixFile, 'ROIs', 'connectivity', ...
            'regionDescriptions', 'weightDescriptions');
        
    end
end

status.reconstruction_network = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_network finished----\n');
