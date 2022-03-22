%% Initialize
testDir = fileparts(mfilename('fullpath'));
assetsDir = fullfile(testDir, 'assets');
testSubjectsDir = fullfile(testDir, 'testSubjects');
toolboxDir = fullfile(fileparts(testDir), 'src');

oldPath = path;
cleanupPath = onCleanup(@() path('',oldPath));

restoredefaultpath;
addpath(genpath(toolboxDir));
addpath(genpath(testDir));

[~,~] = mkdir(testSubjectsDir);

%% Simulate test data
copyfile(fullfile(assetsDir, 'template_nifti.nii.gz'), ...
    fullfile(testSubjectsDir, 'template_nifti.nii.gz'));
testSubjectsFC = simulate_functional_data(testSubjectsDir);
testSubjectsSC = simulate_structural_data(testSubjectsDir);

%% Create and run test suite
import matlab.unittest.TestSuite

import matlab.unittest.parameters.Parameter
import matlab.unittest.selectors.HasParameter;
import matlab.unittest.constraints.StartsWithSubstring
% Option: Add HasParameter('Name','csd') as input variable to
% TestSuite.fromFolder() to select specific tests (where csd is tested).

import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.TestRunner
% Option: runner.addPlugin(StopOnFailuresPlugin)


paramSC = Parameter.fromData('subjectDir', testSubjectsSC);
paramFC = Parameter.fromData('subjectDir', testSubjectsFC);

suiteSCPreprocessing = TestSuite.fromFile( ...
    'structural_pipeline/test_structural_preprocessing.m', ...
    HasParameter('Name', StartsWithSubstring('SC_pre')), ...
    'ExternalParameters', paramSC);

suiteSCReconstructionDiffusion = TestSuite.fromFile( ...
    'structural_pipeline/test_reconstruction_diffusion.m', ...
    HasParameter('Name', 'SC_nonlinearities'), ...    
    'ExternalParameters', paramSC);

suiteSC = TestSuite.fromFolder('structural_pipeline', ...
    HasParameter('Name', 'SC_default'), ...    
    'IncludeSubFolders', true, ...
    'ExternalParameters', paramSC);

suiteFC = TestSuite.fromFolder('functional_pipeline', ...
    'IncludeSubFolders', true, ...
    'ExternalParameters', paramFC);

runner = TestRunner.withTextOutput;   
runner.addPlugin(StopOnFailuresPlugin);

result = runner.run([suiteSCPreprocessing ...
    suiteSCReconstructionDiffusion ...
    suiteSC ...
    suiteFC]);

disp('Results:');
disp(result);

%% Cleanup
rmdir(testSubjectsDir, 's');

assertSuccess(result);