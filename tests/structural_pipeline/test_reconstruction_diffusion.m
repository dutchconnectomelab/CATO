classdef test_reconstruction_diffusion < matlab.unittest.TestCase
    
    properties 
        
        showFig = false;
        
    end
    
    properties (MethodSetupParameter)
        
        subjectDir = {''};
       
    end
    
    properties (TestParameter)
        
        testMethod = {'gqi', 'dti', 'csd'};
        
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
        
        function defaultReconstructionPeaks(testCase, subjectDir, testMethod)
            
            cs = structural_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_SC_ref.json'), ...
                'reconstructionSteps', 'reconstruction_diffusion', ...
                'general.reconstructionMethods', testMethod, ...
                'general.outputDir', 'DWI_processed_test', 'runType', 'overwrite');
            
            diffusionPeaksFile = strrep(cs.reconstruction_diffusion.diffusionPeaksFile, ...
                'METHOD', testMethod);
            obs = load(diffusionPeaksFile);
            ref = load(strrep(diffusionPeaksFile, 'DWI_processed_test', 'CATO_ref'));
     
            % parameters
            numberPeaksTol = 1;
            nansTol = 0.8;
            zerosTol = 0.99;
            angleMedianTol = 3;
            
            % Test all elements are real.
            testCase.verifyTrue(all(isreal(obs.diffusionPeaks)), ...
                'Diffusion peaks contain complex numbers');
            
            % Validate number of reconstructed peaks
            NPeaksTest = sum(any(obs.diffusionPeaks > 0, 2), 3);
            NPeaksRef = sum(any(ref.diffusionPeaks > 0, 2), 3);
            
            testCase.verifyLessThan(std(NPeaksRef - NPeaksTest), numberPeaksTol, ...
                'Number of reconstructed peaks deviates from reference number of peaks');
            % figure; histogram(NPeaksRef - NPeaksTest);
            
            indxNaNTest = any(any(isnan(obs.diffusionPeaks),3), 2);
            indxNaNRef = any(any(isnan(ref.diffusionPeaks),3), 2);
            testCase.verifyGreaterThan(nnz(indxNaNRef == indxNaNTest) ./ length(indxNaNRef), nansTol, ...
                ['Number of diffusion peaks that could not be reconstructed ', ...
                '(NaNs) deviate from reference.']);
            
            indxZerosTest = all(obs.diffusionPeaks(:, :, 1) == 0, 2);
            indxZerosRef = all(ref.diffusionPeaks(:, :, 1) == 0, 2);
            testCase.verifyGreaterThan(nnz(indxZerosRef == indxZerosTest) ./ length(indxZerosRef), zerosTol, ...
                ['Number of diffusion peaks that could not be reconstructed ', ...
                '(zeros) deviate from reference.']);
            
            if contains(testMethod, {'csd', 'gqi'})
                % all voxels
                indx = ~indxZerosTest & ~indxZerosRef;
            elseif contains(testMethod, 'dti')
                % non-crossing-fiber voxels
                indx = sum(any(abs(ref.diffusionPeaks) > 0, 2), 3) == 1;
            else
                error('method not defined');
            end
            
            angDeg = anglePeaks(obs.diffusionPeaks, ref.diffusionPeaks, indx);
            
            testCase.verifyLessThan(median(angDeg), angleMedianTol, ...
                'Median angle between peak deviates from reference');            
            
            if testCase.showFig
                figure; histogram(angDeg);
                xlabel('Angle first peak (in Degrees)');
                s
                figure('Color', 'white');
                hold on;
                for i = 1:size(obs.diffusionPeaks, 3)
                    plotDiffusionPeaks( ...
                        obs.diffusionPeaks(:, :, i), ones(32));
                end
                title([testMethod ' reconstruction'], ...
                    'Interpreter', 'none');
                
            end

        end
        
    end
end