classdef test_reconstruction_functional_network < matlab.unittest.TestCase
    
    properties
    
        reconstructionSteps = {'compute_motion_metrics', 'reconstruction_functional_network'};
    
    end
    
    properties (MethodSetupParameter)
               
        subjectDir = {''};
        
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
        
        function default_reconstruction_functional_network(testCase, subjectDir)
                        
            cf = functional_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_FC_ref.json'), ...
                'reconstructionSteps', testCase.reconstructionSteps, ...
                'general.outputDir', 'fMRI_processed_test');
            
            connectivityMatrixFile = strrep(strrep( ...
                cf.reconstruction_functional_network.connectivityMatrixFile, ...
                'METHOD', 'scrubbed_0.01-0.1'), ...
                'TEMPLATE', 'toyAtlas1');
            
            obs = load(connectivityMatrixFile);
            ref = load(strrep(connectivityMatrixFile, 'fMRI_processed_test', 'CATO_ref'));

            r = corr(squareform(ref.connectivity)', squareform(obs.connectivity)');
            testCase.verifyGreaterThanOrEqual(r, 0.985, 'FC matrix do not correlate enough.')

            
        end
        
        
        
    end
    
end
