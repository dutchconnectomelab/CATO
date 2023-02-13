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
                'general.outputDir', 'fMRI_processed_test', ...
                'runType', 'overwrite');
            
            for reconstructionMethods = {'default_partial', 'default'}

                connectivityMatrixFile = strrep(strrep( ...
                    cf.reconstruction_functional_network.connectivityMatrixFile, ...
                    'METHOD', reconstructionMethods{1}), ...
                    'TEMPLATE', 'toyAtlas1');
                
                obs = load(connectivityMatrixFile);
                ref = load(strrep(connectivityMatrixFile, 'fMRI_processed_test', 'CATO_ref'));
                
                r = corr(squareform(ref.connectivity)', squareform(obs.connectivity)');
                testCase.verifyGreaterThanOrEqual(r, 0.97, ...
                    sprintf('%s FC matrix do not correlate enough.', reconstructionMethods{1}));
                
            end
            
        end
        
        function reconstruction_functional_network_bandpass_filter(testCase, subjectDir)
            % Test the bandpass filter step by adding low and high
            % frequency noise.
            
            % Load original time series.  Voxel in (1,1,1) is a regressor
            % and has normal distributed random signal.
            NT = load_nifti(fullfile(subjectDir, 'CATO_ref/FC_default_fmri.nii.gz'));
            S = NT.vol;
            noiseVoxel = NT.vol(1,1,1,:);                        
            repetitionTime = NT.pixdim(5);
            nScans = size(S,4);
                       
            % Add low and high fequency noise.
            regionPhase = pi .* (linspace(0, 1, 16)'); % 16 brain regions
            regionPhase = regionPhase(randperm(length(regionPhase)));
            noiseLowFreq = 100 * sin([0:(nScans-1)] * 2*pi * (repetitionTime * 0.0066 / 1000)  - regionPhase);
            noiseHighFreq =  100 * sin([0:(nScans-1)] * 2*pi * (repetitionTime * 0.17 / 1000)  - regionPhase);
            
            noise = noiseLowFreq + noiseHighFreq;
            noise = permute(noise, [1 3 4 2]);
            noise = reshape(noise, sqrt(16), sqrt(16), 1, []);
            noise = repelem(noise, 2,2,3);
            noise = noise - min(noise(:));
            
            NT.vol = S + 1*noise;
            NT.vol(1,1,1,:) = noiseVoxel - min(noiseVoxel);
            save_nifti(NT, fullfile(pwd, 'fMRI_processed_test/FC_default_fmri_bandpass_filter.nii.gz'));
                        
            % Run CATO with bandpass filtering ([0.01 - 0.1]).
            cf = functional_pipeline(pwd, ...
                'configurationFile', fullfile(subjectDir, 'config_FC_ref.json'), ...
                'reconstructionSteps', testCase.reconstructionSteps, ...
                'general.outputDir', 'fMRI_processed_test', ...
                'functional_preprocessing.fmriProcessedFile', 'OUTPUTDIR/SUBJECT_fmri_bandpass_filter.nii.gz', ...
                'reconstruction_functional_network.bandpass_filter.filter' , true, ...
                'reconstruction_functional_network_1.bandpass_filter.filter' , true, ...                
                'runType', 'overwrite');
            
            for reconstructionMethods = {'default'}

                connectivityMatrixFile = strrep(strrep( ...
                    cf.reconstruction_functional_network.connectivityMatrixFile, ...
                    'METHOD', reconstructionMethods{1}), ...
                    'TEMPLATE', 'toyAtlas1');           
                
                obs = load(connectivityMatrixFile);
                ref = load(strrep(connectivityMatrixFile, 'fMRI_processed_test', 'CATO_ref'));
                
                r = corr(squareform(ref.connectivity)', squareform(obs.connectivity)');
                testCase.verifyGreaterThanOrEqual(r, 0.98, ...
                    sprintf('%s FC matrices differ too much (in correlation).', reconstructionMethods{1}));
                
                d = mean(abs((squareform(ref.connectivity) - squareform(obs.connectivity))));
                testCase.verifyLessThanOrEqual(d, 0.05, ...
                    sprintf('%s FC matrices differ too much (in absolute distance).', reconstructionMethods{1}));                
                
            end
            
        end        
        
        
        
    end
    
end
