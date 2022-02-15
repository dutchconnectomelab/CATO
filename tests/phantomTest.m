 classdef phantomTest < matlab.unittest.TestCase
    methods(Test)
        function test_phantom_subjects(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tempFixture = testCase.applyFixture(TemporaryFolderFixture);
            
            % load FSL & FreeSurfer paths
            testConfigParams = fileread(fullfile( ...
                fileparts(mfilename('fullpath')), 'test_config.json'));
            testConfigParams = jsondecode(testConfigParams);


            phantomSubjects = testConfigParams.testSubjects;
            for phantom_i = 1:length(phantomSubjects)
                phantomSub = phantomSubjects{phantom_i};
                
                % prepare temporary folder with phantom subject
                subjectDir = fullfile(tempFixture.Folder, phantomSub);
                disp(['The temporary folder: ' subjectDir])
                copyfile(fullfile('tests', 'assets', phantomSub),subjectDir);

                % run CATO
                configFile = fullfile(fileparts(mfilename('fullpath')), ...
                'assets', phantomSub, 'misc/cato_structural_config.json');
                structural_pipeline(subjectDir, ...
                    'configurationFile', configFile, ...
                    'runType', 'overwrite', ...,
                    'general.fslRootDir', testConfigParams.fslRootDir, ...
                    'general.freesurferRootDir', testConfigParams.freesurferRootDir, ...
                    'reconstruction_fibers.NumberOfSeedsPerVoxel', testConfigParams.NumberOfSeedsPerVoxel);

                % check results
                configParams = readConfigFile(configFile);
                methods = configParams.general.reconstructionMethods;
                
                solution = csvread(fullfile(subjectDir, 'misc', 'connectivity.csv'));

                for method_i = 1:length(methods)
                    method = methods{method_i};
                    connectivityFile = [phantomSub  '_connectivity_' ...
                        method '_phantomatlas.mat'];
                    connectivityFile = fullfile(subjectDir, ... 
                        'dwi_processed', connectivityFile);
                    disp(['Checking results in: ', connectivityFile]);
                    res = load(connectivityFile, 'connectivity');
                    connectivity = res.connectivity;
                    c = connectivity(:, :, 3);
                    testCase.assertGreaterThanOrEqual(c, solution/2)
                    testCase.assertGreaterThanOrEqual(1-c, (1 - solution)/2)
                end
            end
        end
    end
end