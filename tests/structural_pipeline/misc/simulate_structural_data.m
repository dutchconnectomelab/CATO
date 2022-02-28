function testSubjectStruct = simulate_structural_data(assetsDir)

oldPwd = pwd;
cleanupPwd = oncleanup(@() cd(oldPwd));

%% Define Nifti template
NT = load_nifti(fullfile(assetsDir, 'template_nifti.nii.gz'));
NT.srow_x = zeros(4,1);
NT.srow_x(1) = 1;
NT.srow_y = zeros(4,1);
NT.srow_y(2) = 1;
NT.srow_z = zeros(4,1);
NT.srow_z(3) = 1;
NT.pixdim = [1 1 1 1 1000 1 1 1]';

%% Default configuration parameters.
testSubject = 'SC_default';

cpDefault = readConfigFile('config_structural_default.json');
cpDefault.general.subject = testSubject;
cpDefault.general.outputDir = 'CATO_ref';
cpDefault.general.templates = 'toyAtlas';
cpDefault.general.templatesDir = 'templateDir';
cpDefault.general.reconstructionMethods =  ...
    {'gqi', 'dti', 'csd', 'gqi_dti', 'csd_dti'};
cpDefault.reconstruction_diffusion.exportNifti.exportNifti = false;

testSubjectStruct = simulate_structural_default(fullfile(assetsDir, testSubject), ...
    cpDefault, NT);

cpDefault = parseConfigParams(cpDefault);

%% SC preprocessing b0 reversed

testSubject = 'SC_pre_b0_reversed';
subjectDir = fullfile(assetsDir, testSubject);
[~,~] = mkdir(subjectDir);
cd(subjectDir);

cpb0rev = readConfigFile('config_structural_default.json');
cpb0rev.general.subject = testSubject;
cpb0rev.general.outputDir = 'CATO_ref';
cpb0rev.general.templates = 'toyAtlas';
cpb0rev.general.templatesDir = 'templateDir';
cpb0rev.structural_preprocessing.preprocessingScript = 'SUBJECT_structural_preprocessing.sh';
cpb0rev.general.reconstructionMethods =  ...
    {'gqi'};
cpb0rev.structural_preprocessing.dwiB0ReversedFile = 'DTI/b0_DTI.nii.gz';

simulate_structural_default(fullfile(assetsDir, testSubject), ...
    cpb0rev, NT);

testSubjectStruct.(testSubject) = subjectDir;

cpb0rev = parseConfigParams(cpb0rev);


dwi = load_nifti(cpb0rev.structural_preprocessing.dwiFile);
nScans = size(dwi.vol,4);
b0 = dwi;
b0.vol = b0.vol(:, :, :, 1);
save_nifti(b0, cpb0rev.structural_preprocessing.dwiB0ReversedFile);

dwi.vol = dwi.vol(:, :, :, [1 1:nScans]); % one b0 scan
save_nifti(dwi, cpb0rev.structural_preprocessing.dwiProcessedFile);

index = [1 2*ones(1,nScans)]';
dlmwrite(cpb0rev.structural_preprocessing.indexFile, index, 'delimiter', ' ');

acqp = [1 0 0 0.085; -1 0 0 0.085];
dlmwrite(cpb0rev.structural_preprocessing.acqpFile, acqp, 'delimiter', ' ');

bvals = [0; dlmread(cpb0rev.structural_preprocessing.rawBvalsFile)];
dlmwrite(cpb0rev.structural_preprocessing.processedBvalsFile, bvals, 'delimiter', ' ');

bvecs = [zeros(1,3); dlmread(cpb0rev.structural_preprocessing.rawBvecsFile)];
dlmwrite(cpb0rev.structural_preprocessing.processedBvecsFile, bvecs, 'delimiter', ' ');

% Check: configParams.structural_preprocessing.b0Scans

fid = fopen(cpb0rev.structural_preprocessing.preprocessingScript, 'w+');
fprintf(fid, 'echo "Example preprocessing script"');
fclose(fid);
system(['chmod 700 ' cpb0rev.structural_preprocessing.preprocessingScript]);


%% SC preprocessing b0 reversed

testSubject = 'SC_pre_reversed';
subjectDir = fullfile(assetsDir, testSubject);
[~,~] = mkdir(subjectDir);
cd(subjectDir);

cpRev = readConfigFile('config_structural_default.json');
cpRev.general.subject = testSubject;
cpRev.general.outputDir = 'CATO_ref';
cpRev.general.templates = 'toyAtlas';
cpRev.general.templatesDir = 'templateDir';
cpRev.structural_preprocessing.preprocessingScript = 'SUBJECT_structural_preprocessing.sh';
cpRev.general.reconstructionMethods =  ...
    {'gqi'};
cpRev.structural_preprocessing.dwiReversedFile = 'DTI/reversed_DTI.nii.gz';

simulate_structural_default(subjectDir, ...
    cpRev, NT);

testSubjectStruct.(testSubject) = subjectDir;

cpRev = parseConfigParams(cpRev);

dwi = load_nifti(cpRev.structural_preprocessing.dwiFile);
nScans = size(dwi.vol,4);
save_nifti(dwi, cpRev.structural_preprocessing.dwiReversedFile);

dwi.vol = cat(4, dwi.vol, dwi.vol);
save_nifti(dwi, cpRev.structural_preprocessing.dwiProcessedFile);

index = [ones(1,nScans) 2*ones(1,nScans)]';
dlmwrite(cpRev.structural_preprocessing.indexFile, index, 'delimiter', ' ');

acqp = [1 0 0 0.085; -1 0 0 0.085];
dlmwrite(cpRev.structural_preprocessing.acqpFile, acqp, 'delimiter', ' ');

bvals = dlmread(cpRev.structural_preprocessing.rawBvalsFile);
bvals = [bvals; bvals];
dlmwrite(cpRev.structural_preprocessing.processedBvalsFile, bvals, 'delimiter', ' ');

bvecs = dlmread(cpRev.structural_preprocessing.rawBvecsFile);
bvecs = [bvecs; bvecs];
dlmwrite(cpRev.structural_preprocessing.processedBvecsFile, bvecs, 'delimiter', ' ');

% Check: configParams.structural_preprocessing.b0Scans
fid = fopen(cpRev.structural_preprocessing.preprocessingScript, 'w+');
fprintf(fid, 'echo "Example preprocessing script"');
fclose(fid);
system(['chmod 700 ' cpRev.structural_preprocessing.preprocessingScript]);


%% SC_nonlinearities

testSubject = 'SC_nonlinearities';
subjectDir = fullfile(assetsDir, testSubject);
[~,~] = mkdir(subjectDir);
cd(subjectDir);

cpNonLin = readConfigFile('config_structural_default.json');
cpNonLin.general.subject = testSubject;
cpNonLin.general.outputDir = 'CATO_ref';
cpNonLin.general.reconstructionMethods =  ...
    {'gqi', 'dti', 'csd'};
cpNonLin.reconstruction_diffusion...
    .gradientNonlinearities.correctNonlinearities = true;
cpNonLin.reconstruction_diffusion...
    .gradientNonlinearities.nonlinearitiesFile = 'nonlinearities.nii.gz';
cpNonLin.general.templates = 'toyAtlas';
cpNonLin.general.templatesDir = 'templateDir';
saveConfigFile('config_SC_ref.json', cpNonLin);

subjectDir = fullfile(assetsDir, testSubject);

simulate_structural_default(subjectDir, ...
    cpNonLin, NT);

cpNonLin = parseConfigParams(cpNonLin);

testSubjectStruct.(testSubject) = subjectDir;

% create non-linearities file (rotate xy-plane)
nonlinearities = [-1 -1 0; 1 -1 0; 0 0 0];
nonlinearities = nonlinearities(:);

NL = NT;
NL.vol = repmat(permute(nonlinearities, [2 3 4 1]), ...
    [size(NL.vol, 1), size(NL.vol, 2), size(NL.vol, 3) 1]);

save_nifti(NL, fullfile(subjectDir, cpNonLin.reconstruction_diffusion...
    .gradientNonlinearities.nonlinearitiesFile));

% Create rotated diffusionPeaksFile
for im = 1:length(cpNonLin.general.reconstructionMethods)
    m = cpNonLin.general.reconstructionMethods{im};
    
    ref = load(fullfile(subjectDir, strrep(cpNonLin...
        .reconstruction_diffusion.diffusionPeaksFile, ...
        'METHOD', m)));
    diffusionPeaksNew = nan(size(ref.diffusionPeaks));
    for i = 1:size(ref.diffusionPeaks, 3)
        diffusionPeaksNew(:, :, i) = ...
            ref.diffusionPeaks(:, :, i) * (eye(3) + ...
            reshape(nonlinearities, [3 3]));
    end
    ref.diffusionPeaks = diffusionPeaksNew;
    save(fullfile(subjectDir, strrep(cpNonLin... 
        .reconstruction_diffusion.diffusionPeaksFile, ...
        'METHOD', m)), '-struct', 'ref');
    
end
 
