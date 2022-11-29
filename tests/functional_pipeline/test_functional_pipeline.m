classdef test_functional_pipeline < matlab.unittest.TestCase
    
    properties
        reconstructionSteps = {'compute_motion_metrics'};
    end
    
    properties (MethodSetupParameter)
        
        subjectDir = {''};
       
    end
    
    
    methods (TestClassSetup)
        
        
    end
    
    methods (TestMethodSetup)
        
        function prepare_data(testCase,subjectDir)
            
            import matlab.unittest.fixtures.WorkingFolderFixture;
            testCase.applyFixture(WorkingFolderFixture);
            
            [~, subjectName] = fileparts(subjectDir);
            subjectDirTMP = fullfile(pwd, subjectName);
            copyfile(subjectDir, subjectDirTMP);
            copyfile(fullfile(subjectName, 'CATO_ref'), ...
                fullfile(subjectName, 'fMRI_processed_test'));
            
            import matlab.unittest.fixtures.CurrentFolderFixture
            testCase.applyFixture(CurrentFolderFixture ...
                (subjectDirTMP));
            
        end
        
    end
    
    methods(Test)
        
        function errorHandling(testCase, subjectDir)
            
            fid = fopen('test.par', 'w+');
            fclose(fid);
                        
            cf = functional_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_FC_ref.json'), ...
                'reconstructionSteps', testCase.reconstructionSteps, ...
                'functional_preprocessing.motionParametersFile', 'test.par', ...
                'general.outputDir', 'fMRI_processed_test');            
            
            testCase.assertTrue(isfile(cf.general.logFile));
            logText = textread(cf.general.logFile, '%s', 'delimiter', '####');
            
            testCase.verifyTrue(any(cellfun(@(x) contains(x, 'Error in: compute_motion_metrics (Line'), logText)));
            
            testCase.assertTrue(isfile(strrep(cf.general.statusFile, 'STATUS', 'error')));

            status = readConfigFile(strrep(cf.general.statusFile, 'STATUS', 'error'));
            
            testCase.verifyTrue(isequal(status.general, 'error'));
            testCase.verifyTrue(isequal(status.compute_motion_metrics, 'error'));
            
                  
        end
        
    end
end