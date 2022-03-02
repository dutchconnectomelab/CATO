classdef test_structural_pipeline < matlab.unittest.TestCase
    
    
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
                fullfile(subjectName, 'DWI_processed_test'));
            
            import matlab.unittest.fixtures.CurrentFolderFixture
            testCase.applyFixture(CurrentFolderFixture ...
                (subjectDirTMP));
            
        end
        
    end
    
    methods(Test)
        
        function errorHandling(testCase, subjectDir)
            
            fid = fopen('test.nii.gz', 'w+');
            fclose(fid);
            
            cs = structural_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_SC_ref.json'), ...
                'reconstructionSteps', 'structural_preprocessing', ...
                'structural_preprocessing.preprocessingScript', fullfile(subjectDir, 'SUBJECT_structural_preprocessing.sh'), ...
                'general.outputDir', 'DWI_processed_test', ...
                'structural_preprocessing.dwiFile', 'test.nii.gz', ...
                'runType', 'overwrite');
            
            testCase.assertTrue(isfile(cs.general.logFile));
            logText = textread(cs.general.logFile, '%s', 'delimiter', '####');
            
            testCase.verifyTrue(any(cellfun(@(x) contains(x, 'Error in: structural_preprocessing (Line'), logText)));
            
            testCase.assertTrue(isfile(strrep(cs.general.statusFile, 'STATUS', 'error')));

            status = readConfigFile(strrep(cs.general.statusFile, 'STATUS', 'error'));
            
            testCase.verifyTrue(isequal(status.general, 'error'));
            testCase.verifyTrue(isequal(status.structural_preprocessing, 'error'));
            
                  
        end
        
    end
end