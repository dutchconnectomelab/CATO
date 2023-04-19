classdef test_reconstruction_fibers < matlab.unittest.TestCase
    
    
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
                fullfile(subjectName, 'DWI_processed_test'));
            
            import matlab.unittest.fixtures.CurrentFolderFixture
            testCase.applyFixture(CurrentFolderFixture ...
                (subjectDirTMP));
            
        end
        
    end
    
    methods(Test)
        
        function testOrientation(testCase, subjectDir)
            
            [~, subjectName] = fileparts(subjectDir);
            cp = readConfigFile(fullfile(subjectDir, 'config_SC_ref.json'));
            cp.general.subject = subjectName;
            cp = parseConfigParams(cp);
            
            sform = eye(3);
            sform(:, :, 2) = [0 -1 0; 1 0 0; 0 0 1];
            sform(:, :, 3) = ...
                            [[   -2.0000    0.0000   -0.0000]; ...
                            [    0.0000    1.7560   -0.9574]; ...
                            [         0    0.9574    1.7560]];
            sform(:, :, 4) = [0 0 -1; 0 -1 0; -1 0 0];
            ref = {'RAS', 'ALS', 'LAS', 'IPL'}; 
            
            for iPerm = 1:size(sform,3)
                
                NT = load_nifti(cp.structural_preprocessing.segmentationFile);
                NT.srow = sform(:, :,iPerm);
                NT.srow_x = [sform(1,:, iPerm) 0]';
                NT.srow_y  = [sform(2,:, iPerm) 0]';
                NT.srow_z  = [sform(3,:, iPerm) 0]';
                NT.pixdim = [1 1 1 1 1000 1 1 1]';
                save_nifti(NT, strrep(cp.structural_preprocessing.segmentationFile, ...
                    'CATO_ref', 'DWI_processed_test'));
                    
                cs = structural_pipeline(pwd, ...
                    'configurationFile', fullfile(subjectDir, 'config_SC_ref.json'), ...
                    'general.outputDir', 'DWI_processed_test', ...
                    'reconstructionSteps', 'reconstruction_fibers', ...
                    'general.reconstructionMethods', 'csd', ...
                    'reconstruction_fibers.NumberOfSeedsPerVoxel', 1, 'runtype', 'overwrite');                
                
                [~, header] = readTrk(strrep(cs.reconstruction_fibers.fiberFile, ...
                    'METHOD', 'csd'));
                
                obs = strtrim(strip(header.voxel_order, char(0)));
                
                fprintf('Reference orientation: %s\n', ...
                    ref{iPerm});

                % Show the orientation according to mri_info
                fprintf('mri_info:');
                system(sprintf(['source ~/.bashrc; ' ...
                    'mri_info --orientation %s | tail -n1'], ...
                    strrep(cp.structural_preprocessing.segmentationFile, ...
                    'CATO_ref', 'DWI_processed_test')));
                
                testCase.verifyEqual(obs, ref{iPerm}, 'Orientation not correct.')
                
            end
            
        end
        
    end
end