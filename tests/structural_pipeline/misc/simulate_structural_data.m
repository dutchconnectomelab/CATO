function testSubjectStruct = create_simulated_data(assetsDir)

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

%% SC_nonlinearities

testSubject = 'SC_nonlinearities';

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
 
