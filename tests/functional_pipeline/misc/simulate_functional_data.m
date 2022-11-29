function testSubjectStruct = simulate_functional_data(assetsDir)
% SIMULATE_FUNCTIONAL_DATA simulate rs-fMRI data
%

testSubject = 'FC_default';
testSubjectDir = fullfile(assetsDir, testSubject);
[~, ~] = mkdir(testSubjectDir);

testSubjectStruct = struct(testSubject, testSubjectDir);

oldPath = pwd;
cleanupPath = onCleanup(@() cd(oldPath));

cd(fullfile(assetsDir, testSubject));

% PARAMETERS
% Simulated signal has dimensions:
% sqrt(nRegions) x sqrt(nRegions) x 3 x nScans
nRegions = 16;
nScans = 500;
repetitionTime = 750;


% Load default configuration parameters.
configParams = readConfigFile('config_functional_default.json');
configParams.general.subject = testSubject;
configParams.general.outputDir = 'CATO_ref';
configParams.general.templates = {'toyAtlas1'};
configParams.general.templatesDir = 'templateDir';

configParams.reconstruction_functional_network.methodDescription = 'default';
configParams.reconstruction_functional_network.bandpass_filter.filter = false;
configParams.reconstruction_functional_network.scrubbing.scrubbing = false;
configParams.reconstruction_functional_network.regression.regressionMask = 999;

configParams.reconstruction_functional_network_1 = ...
    configParams.reconstruction_functional_network;
configParams.reconstruction_functional_network_1.methodDescription = 'default_partial';
configParams.reconstruction_functional_network_1.bandpass_filter.filter = false;
configParams.reconstruction_functional_network_1.scrubbing.scrubbing = false;
configParams.reconstruction_functional_network_1.regression.regressionMask = 999;
configParams.reconstruction_functional_network_1.reconstructionMethod = 'partialcorr';

saveConfigFile('config_FC_ref.json', configParams);

configParams = parseConfigParams(configParams);

% Setup default Nifti template:
NT = load_nifti(fullfile(assetsDir, 'template_nifti.nii.gz'));
NT.srow_x = zeros(4,1);
NT.srow_x(1) = 1;
NT.srow_y = zeros(4,1);
NT.srow_y(2) = 1;
NT.srow_z = zeros(4,1);
NT.srow_z(3) = 1;
NT.pixdim = [1 1 1 1 repetitionTime 1 1 1]';

% INITIALIZATION
[~, ~] = mkdir(configParams.general.outputDir);

%% Simulate rs-fMRI data

% Make covariance matrix and model rs-fMRI data as a multivariate normal
% distribution.
refCov = rand(nRegions,nRegions);
refCov(end+1:end+7, end+1:end+7) = eye(7);
refCov = normalize(refCov, 'norm');
refCov = refCov' * refCov;

S = 1000*mvnrnd(10 * ones(nRegions+7,1), refCov,nScans);

noise = S(:,nRegions+1);
rotationParams = S(:,nRegions+2:nRegions+4);
translationParams = S(:,nRegions+5:nRegions+7);
S = S(:, 1:nRegions);

% partial correlations:
% Sp = zeros(size(S));
% X = [rotationParams translationParams];
% X = [X [zeros(1,size(X,2)); diff(X', 1, 2)'] noise];
% X = bsxfun(@rdivide, X', mean(abs(X'), 2))';
% for i = 1:nRegions
%     SpMdl = fitlm(X, S(:,i), 'intercept', false);
%     Sp(:,i) = SpMdl.Residuals.Raw - mean(SpMdl.Residuals.Raw);
% end

data_partial.connectivity = partialcorr(S);
data_partial.connectivity(eye(nRegions)>0) = 0;

S = permute(S, [2 3 4 1]);
S = reshape(S, sqrt(nRegions), sqrt(nRegions), 1, []);
S = repelem(S, 2,2,3);

% Create one voxel with noise to regress out in regression step.
S(1,1,1,:) = noise;

NT.vol = S;
save_nifti(NT, configParams.functional_preprocessing.fmriProcessedFile);

data.connectivity = refCov(1:nRegions, 1:nRegions);
data.connectivity(eye(nRegions)>0) = 0;

save(strrep(strrep(configParams.reconstruction_functional_network.connectivityMatrixFile, ...
    'METHOD', 'default'), ...
    'TEMPLATE', 'toyAtlas1'), ...
    '-struct', 'data');

save(strrep(strrep(configParams.reconstruction_functional_network_1.connectivityMatrixFile, ...
    'METHOD', 'default_partial'), ...
    'TEMPLATE', 'toyAtlas1'), ...
    '-struct', 'data_partial');

%% Create segmentationFile

segmentation = 1:nRegions;
segmentation = reshape(segmentation, [sqrt(nRegions), sqrt(nRegions)]);
segmentation = repelem(segmentation, 2,2,3);
segmentation(1) = 999;

NT.vol = segmentation;

save_nifti(NT, configParams.functional_preprocessing.segmentationFile);
save_nifti(NT, strrep(configParams.parcellation.parcellationFile, ...
    'TEMPLATE', 'toyAtlas1'));
save_nifti(NT, strrep(configParams.parcellation.parcellationFile, ...
    'TEMPLATE', 'toyAtlas2'));

%% Create motionParametersFile

dlmwrite(configParams.functional_preprocessing.motionParametersFile, ...
    [rotationParams translationParams]);

%% Create template directories

templates = {'toyAtlas1', 'toyAtlas2'};

for iT = 1:2

ROIs = (1:nRegions)';
names = strcat('region_', string(num2cell(1:nRegions)))';

x = table(ROIs, names, ...
    rand(nRegions,1), rand(nRegions,1), ...
    rand(nRegions,1), rand(nRegions,1));

[~, ~] = mkdir('templateDir');
[~, ~] = mkdir(fullfile('templateDir', templates{iT}));

writetable(x, strrep('templateDir/TEMPLATE/TEMPLATE.annot.ctab', ...
    'TEMPLATE', templates{iT}), 'delimiter', '\t', ...
    'WriteVariableNames', false, ...
    'QuoteStrings', false, ...
    'fileType', 'text');

dlmwrite(strrep('templateDir/TEMPLATE/ROIs_TEMPLATE.txt', ...
    'TEMPLATE', templates{iT}), ROIs);

end