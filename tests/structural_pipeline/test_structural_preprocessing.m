classdef test_structural_preprocessing < matlab.unittest.TestCase
    
    
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
        
        function defaultPreprocessing(testCase, subjectDir)
            
            cs = structural_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_SC_ref.json'), ...
                'reconstructionSteps', 'structural_preprocessing', ...
                'structural_preprocessing.preprocessingScript', fullfile(subjectDir, 'SUBJECT_structural_preprocessing.sh'), ...
                'general.outputDir', 'DWI_processed_test', 'runType', 'overwrite');
            
            indexFile = cs.structural_preprocessing.indexFile;
            obs = dlmread(indexFile);
            ref = dlmread(strrep(indexFile, 'DWI_processed_test', 'CATO_ref'));
            
            testCase.verifyEqual(obs, ref);            
            
            acqpFile = cs.structural_preprocessing.acqpFile;
            obs = dlmread(acqpFile);
            ref = dlmread(strrep(acqpFile, 'DWI_processed_test', 'CATO_ref'));            
            
            testCase.verifyEqual(obs, ref);
            
            processedBvalsFile = cs.structural_preprocessing.processedBvalsFile;
            obs = dlmread(processedBvalsFile);
            ref = dlmread(strrep(processedBvalsFile, 'DWI_processed_test', 'CATO_ref'));            
            
            testCase.verifyEqual(obs, ref);
            
            processedBvecsFile = cs.structural_preprocessing.processedBvecsFile;
            obs = dlmread(processedBvecsFile);
            ref = dlmread(strrep(processedBvecsFile, 'DWI_processed_test', 'CATO_ref'));            
            
            testCase.verifyEqual(obs, ref);
            
            dwiProcessedFile = cs.structural_preprocessing.dwiProcessedFile;
            obs = load_nifti(dwiProcessedFile);
            ref = load_nifti(strrep(dwiProcessedFile, 'DWI_processed_test', 'CATO_ref'));            
            
            testCase.verifyEqual(obs, ref);            
     
        

        end
        
    end
end