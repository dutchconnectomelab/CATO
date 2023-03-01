function commit_filter(configParams)
% COMMIT_FILTER   Filter the fiber cloud using COMMIT2
%

%% Initialization
status.commit_filter = 'running';
updateStatus(configParams.general.statusFile, status);
fprintf('---commit_filter started----\n');

dwiProcessedFile = configParams.structural_preprocessing.dwiProcessedFile;

fiberFile = configParams.reconstruction_fibers.fiberFile;
fiberPropertiesFile = configParams.reconstruction_fiber_properties.fiberPropertiesFile;

regionPropertiesFile = configParams.collect_region_properties.regionPropertiesFile;

processedBvalsFile = configParams.structural_preprocessing.processedBvalsFile;
processedBvecsFile = configParams.structural_preprocessing.processedBvecsFile;
bValueScalingTol = configParams.general.bValueScalingTol;
bValueZeroThreshold = configParams.general.bValueZeroThreshold;


reconstructionMethods = lower(configParams.commit_filter.reconstructionMethods);
templates = configParams.commit_filter.templates; % Run COMMIT-filter only on a selection of templates
filteredFiberFile = configParams.commit_filter.filteredFiberFile;
filteredFiberPropertiesFile = configParams.commit_filter.filteredFiberPropertiesFile;
schemeFile = configParams.commit_filter.schemeFile;
intermediateConnectomeFile = configParams.commit_filter.intermediateConnectomeFile;
pythonInterpreter = configParams.commit_filter.pythonInterpreter;
commitScriptFile = configParams.commit_filter.commitScriptFile;
fiberWeightsFile = configParams.commit_filter.fiberWeightsFile;
filteredConnectivityMatrixFile = configParams.commit_filter.filteredConnectivityMatrixFile;
outputCommitDir = configParams.commit_filter.outputCommitDir;
lambda = configParams.commit_filter.lambda;
subjectDir = pwd;


%% Convert bvals and bvecs to single scheme file

% DWI.scheme has the following organization:
% VERSION: BVECTOR
% -0.000000	-0.000000	1.000000	0.000000
% 0.971803	-0.194807	-0.132847	700.000000

gtab = load_gtab(processedBvalsFile, processedBvecsFile, bValueZeroThreshold, bValueScalingTol);


[~, ~] = mkdir(fileparts(schemeFile));
fid = fopen(schemeFile, 'w+');

if fid == -1
     error('CATO:commit_filter:cannotCreateFile', 'Cannot create ''%s''.', schemeFile);
end

fprintf(fid, 'VERSION: BVECTOR\n');
fprintf(fid, '%.6f\t%.6f\t%.6f\t%.6f\n', [gtab.bvecs gtab.bvals]');

fclose(fid);


%% Reconstruct networks for each template and each method.
for iTemplate = 1:length(templates)
    thisTemplate = templates{iTemplate};
    
    for iMethod = 1:length(reconstructionMethods)
        thisMethod = reconstructionMethods{iMethod};
        
        thisFiberFile =  strrep(fiberFile, ...
            'METHOD', thisMethod);
        thisFilteredFiberFile =  strrep(filteredFiberFile, ...
            'METHOD', thisMethod);
        thisFiberPropertiesFile =  strrep(strrep(fiberPropertiesFile, ...
            'METHOD', thisMethod), ...
            'TEMPLATE', thisTemplate);
        thisFilteredFiberPropertiesFile =  strrep(strrep(filteredFiberPropertiesFile, ...
            'METHOD', thisMethod), ...
            'TEMPLATE', thisTemplate);        
        thisRegionPropertiesFile = strrep(regionPropertiesFile, ...
            'TEMPLATE', thisTemplate);        

        fprintf('%s - %s\n', thisTemplate, upper(thisMethod));

        % Open reconstructed fiber file.
        % This script assumes that the fibercloud data is small enough to
        % fit into memory, which is quite a liberal assumption.
        fprintf('Load fibers...');
        [fibers, headerFiberFile] = readTrk(thisFiberFile);
        fprintf(' done\n');

        % Fiber properties
        data = load(thisFiberPropertiesFile);
        fiberProperties = data.fiberProperties;
        propertyDescriptions = data.propertyDescriptions;
        
        % Filter fibers that are in connectome and group fibers into
        % connection bundles. 

        % To keep things simple we want each fiber only one time in the
        % fiber cloud file. We will order fibers by the longest connection
        % they contribute to.

        fprintf('Filter fibers...');
        % Sort such that first occurances of a fiber are the longest
        [~, indx] = sort(fiberProperties(:,6), 'descend');
        fiberProperties = fiberProperties(indx,:);

        % Filter the first occurances of a fiber
        [~, indx] = unique(fiberProperties(:,1), 'stable');
        fiberProperties = fiberProperties(indx,:);

        % Sort/group the filtered fibers by connection
        [~, indx] = sort(fiberProperties(:,2), 'ascend');
        fibers = fibers(indx);
        fprintf(' done\n');

        fprintf('Write intermediate files for Python COMMIT script...');
        voxelSize = headerFiberFile.voxel_size;
        headerFiberFile.voxel_order = strtrim(headerFiberFile.voxel_order);
        writeFibers(fibers, thisFilteredFiberFile, voxelSize, headerFiberFile);

        % Create connectome file for COMMIT script with the updated NOS.
        data = load(thisRegionPropertiesFile);
        regionDescriptions = data.regionDescriptions;        
        nRegions = length(regionDescriptions);
        
        assert(max(fiberProperties(:, 2)) <= 0.5*nRegions*(nRegions-1), ...
            ['The number of regions does not match between ', ...
            'fiberPropertiesFile (%s) and regionPropertiesFile (%s).'], ...
            thisFiberPropertiesFile, thisRegionPropertiesFile);

        W = accumarray(fiberProperties(:, 2), ones(size(fiberProperties,1),1), ...
            [0.5*nRegions*(nRegions-1) 1]);
        W = squareform(W);        
        dlmwrite(intermediateConnectomeFile, W, ',');
        fprintf(' done\n');

        % Run Python script to perform COMMIT and COMMIT2
        cmd = [pythonInterpreter ' "' commitScriptFile '"' ...
            ' --dwiProcessedFile=' dwiProcessedFile, ...
            ' --dwiSchemeFile=' schemeFile, ...
            ' --subjectDir=' subjectDir, ...
            ' --fiberFile=' thisFilteredFiberFile, ...
            ' --outputCommitDir=' outputCommitDir, ...
            ' --regLambda=' num2str(lambda), ...
            ' --intermediateConnectomeFile=' intermediateConnectomeFile];

        fprintf('COMMIT script:\n%s\n', commitScriptFile);
        fprintf('Run COMMIT script...');
        exitCode = system(cmd);

        if exitCode ~=0
            error('CATO:commit_filter:errorInCOMMITScript', ...
                'Error in executing COMMIT script: \n %s', ...
                cmd);
        end
        fprintf(' done\n');

        fprintf('Apply COMMIT filter to connectome matrices...\n');
        % Filter with streamline weights
        commitWeights = dlmread(fiberWeightsFile);
        
        % Create updated fiber properties file
        fiberProperties = fiberProperties(commitWeights > 0,:);
        save(thisFilteredFiberPropertiesFile, 'fiberProperties', 'propertyDescriptions');
        
        % Run network reconstruction step
        configParams.reconstruction_fiber_properties.fiberPropertiesFile = filteredFiberPropertiesFile;
        configParams.general.templates = {thisTemplate};
        configParams.general.reconstructionMethods = {thisMethod};
        configParams.reconstruction_network.connectivityMatrixFile = filteredConnectivityMatrixFile;
        reconstruction_network(configParams);
        
        fprintf('done.\n')
    end 

end

%% Clean up

status.commit_filter = 'finished';
updateStatus(configParams.general.statusFile, status);

fprintf('---commit_filter finished----\n');


