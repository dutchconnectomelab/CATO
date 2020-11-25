function reconstruction_fiber_properties(configParams)
% RECONSTRUCTION_FIBER_PROPERTIES   Find fiber segments and calculate
%                                   associated segment properties.
%
%   reconstruction_fiber_properties(CONFIGPARAMS) find all fiber segments
%   that connection ROIs and calculate the associated segment properties
%   based on the parameters specified in CONFIGPARAMS.

%% Initialize
status.reconstruction_fiber_properties = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_fiber_properties started----\n');

% Variables
diffusionMeasuresFile = configParams.reconstruction_diffusion.diffusionMeasuresFile;
parcellationFile = configParams.parcellation.parcellationFile;
ROIsFile = configParams.general.ROIsFile;
fiberFile = configParams.reconstruction_fibers.fiberFile;
fiberPropertiesFile = configParams.reconstruction_fiber_properties.fiberPropertiesFile;
includeGMVoxelsFlag = configParams.reconstruction_fiber_properties.includeGMVoxelsFlag;
reconMethods = lower(configParams.general.reconstructionMethods);
templates = configParams.general.templates;

% Load data
data = load(diffusionMeasuresFile);
diffusionMeasures = data.diffusionMeasures;
diffusionMeasures = reshape(diffusionMeasures, [], size(diffusionMeasures, 4));
weightDescriptions = data.weightDescriptions;

%% Compute fiber properties

for iTemplate = 1:length(templates)
for iMethod = 1:length(reconMethods) 
    thisMethod = reconMethods{iMethod};
    thisTemplate = templates{iTemplate};

    % Prepare dynamic variables
    thisParcellationFile = strrep(parcellationFile, 'TEMPLATE', thisTemplate);
    thisROIsFile = strrep(ROIsFile, 'TEMPLATE', thisTemplate);
    thisFiberFile = strrep(fiberFile, 'METHOD', thisMethod);
    thisfiberPropertiesFile =  strrep(strrep(fiberPropertiesFile, ...
                               'METHOD', thisMethod), ...
                               'TEMPLATE', thisTemplate);

    % Load dynamic data
    parcellation = load_nifti(thisParcellationFile);
    ROIsList = dlmread(thisROIsFile);

    % Compute fiber properties
    [fiberProperties, propertyDescriptions] = getFiberPropertiesFromFile(thisFiberFile, ...  
        ROIsList, parcellation, diffusionMeasures, weightDescriptions, includeGMVoxelsFlag); 
    
    % Save properties
    save(thisfiberPropertiesFile, 'fiberProperties', 'propertyDescriptions');

end
end

status.reconstruction_fiber_properties = 'finished';
updateStatus(configParams.general.statusFile, status);
fprintf('---reconstruction_fiber_properties finished----\n');




