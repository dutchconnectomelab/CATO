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

import matlab.unittest.parameters.Parameter
import matlab.unittest.selectors.HasParameter;
% Option: Add HasParameter('Name','csd') as input variable to to
% TestSuite.fromFolder() to select specific tests (where csd is tested).

import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.TestRunner
% Option: runner.addPlugin(StopOnFailuresPlugin)


paramSC = Parameter.fromData('subjectDir', testSubjectsSC);
paramFC = Parameter.fromData('subjectDir', testSubjectsFC);

import matlab.unittest.TestSuite
suiteSC = TestSuite.fromFolder('structural_pipeline', ...
    'IncludeSubFolders', true, ...
    'ExternalParameters', paramSC);

suiteFC = TestSuite.fromFolder('functional_pipeline', ...
    'IncludeSubFolders', true, ...
    'ExternalParameters', paramFC);

runner = TestRunner.withTextOutput;   
runner.addPlugin(StopOnFailuresPlugin);
result = runner.run([suiteSC suiteFC])

%% Cleanup

rmdir(testSubjectsDir, 's');
