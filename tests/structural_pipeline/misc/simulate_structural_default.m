function testSubjectStruct = simulate_structural_default(testSubjectDir, cpDefault, NT)

S0 = 1500;
noise_level = 0;
showFig = false;

gtab.bvecs = [0 0 0; 0.1639 0.5115 0.8435;0.1176 -0.5388 0.8342; ...
    0.5554 0.8278 -0.0797;-0.4804 0.8719 0.0948; ...
    0.9251 -0.0442 0.3772;0.7512 -0.0273 -0.6596; ...
    0.1655 -0.0161 0.9861;0.6129 -0.3427 0.7120; ...
    0.6401 0.2747 0.7175;-0.3724 -0.3007 0.8780; ...
    -0.3451 0.3167 0.8835;0.4228 0.7872 0.4489; ...
    0.0441 0.9990 0.0089;-0.1860 0.8131 0.5515; ...
    0.8702 0.4606 0.1748;-0.7239 0.5285 0.4434; ...
    -0.2574 -0.8032 0.5372;0.3515 -0.8292 0.4346; ...
    -0.7680 -0.4705 0.4346;0.8261 -0.5384 0.1660; ...
    0.9852 -0.0420 -0.1660];

gtab.bvals = [0; 1500*ones(21,1)];

testSubjectStruct = struct(cpDefault.general.subject, testSubjectDir);

[~, ~] = mkdir(testSubjectDir);

oldPath = pwd;
cleanupPath = onCleanup(@() cd(oldPath));
cd(testSubjectDir);

[~,~] = mkdir(cpDefault.general.outputDir);

saveConfigFile('config_SC_ref.json', cpDefault);
cpDefault = parseConfigParams(cpDefault);

%% Simulate DWI signal

%The next lines implement the simulation of the DW-MRI signal
S=zeros(32,32,size(gtab.bvecs,1));
fO1 = zeros(32,32,3);
fO2 = zeros(32,32,3);
for i=1:size(gtab.bvecs,1)
    fprintf(1,'.');
    for iR=1:32
        for y=1:32
            
            f1_flag=0;
            f2_flag=0;
            
            %Define the underlying fiber orientation of fiber bundle #1
            if iR*iR+y*y>16*16 & iR*iR+y*y<32*32
                v=[y/iR -1 0];v=v/sqrt(v*v');
                fiber_orientation1=v;
                f1_flag=1;
            end
            
            %Define the underlying fiber orientation of fiber bundle #2 (i.e. we have 2-fiber crossing)
            if iR<y+10 & iR>y-10
                fiber_orientation2=[sqrt(2)/2 sqrt(2)/2 0];
                f2_flag=1;
            end
            
            %Here we fill the rest of the voxels
            if f1_flag==0 & f2_flag==1
                fiber_orientation1=fiber_orientation2;
            elseif f1_flag==1 & f2_flag==0
                fiber_orientation2=fiber_orientation1;
            elseif f1_flag==0 & f2_flag==0
                fiber_orientation1 = [0 0 0]; % sqrt(3)/3*ones(1,3); % [0 0 1];
                fiber_orientation2 = [0 0 0]; % sqrt(3)/3*ones(1,3); % [0 0 1];
            end
            
            fO1(iR, y, :) = fiber_orientation1;
            fO2(iR, y, :) = fiber_orientation2;
            
            S(iR,y,i)=S0 * ...
                (SimulateDWMRI(fiber_orientation1, gtab.bvecs(i, :))+ ...
                SimulateDWMRI(fiber_orientation2, gtab.bvecs(i, :))) / 2;
        end
    end
end

S(:, :, gtab.bvals<1) = S0;
S = permute(S, [1 2 4 3]);
S = cat(3, S, S, S);

fO1v = reshape(fO1, 32*32, []);
fO2v = reshape(fO2, 32*32, []);

if showFig
    figure('Color', 'white');
    hold on;
    V = ones(32,32);
    plotDiffusionPeaks(fO1v, V);
    plotDiffusionPeaks(fO2v, V);
    title('ground truth');
end

% Simulate noise
Sv = reshape(S, 32*32*3, []);

noise1 = noise_level*std(Sv).*randn(size(Sv));
noise2 = noise_level*std(Sv).*randn(size(Sv));
Sv = sqrt((Sv + noise1).^2 + noise2.^2);

S = reshape(Sv, 32, 32, 3, []);



if showFig
    % Show an example diffusion peak
    % load peak_indices reconstruction basis    
    
    indxVoxel = 17*32+17;
    P = [gtab.bvecs(2:end,:); -gtab.bvecs(2:end,:)];
    
    figure('color', 'white');
    hold on;
    k = boundary(P,0);
    p = [Sv(indxVoxel,2:end)'; Sv(indxVoxel,2:end)'];
    p = -log(p / S0);
    Pp = P .* p;
    trisurf(k,Pp(:,1),Pp(:,2),Pp(:,3),p, 'FaceColor', 'interp', ...
        'EdgeColor','none','SpecularStrength',0);
    axis square
    daspect([1 1 1])
    view(3); axis tight
    
    pm = max(Pp(:))*1.5;
    h = line(pm.*[-fO1v(indxVoxel,1) fO1v(indxVoxel,1)], ...
        pm.*[-fO1v(indxVoxel,2) fO1v(indxVoxel,2)], ...
        pm.*[-fO1v(indxVoxel,3) fO1v(indxVoxel,3)], ...
        'Color', [0 0 0], 'LineWidth', 2);
    h = line(pm.*[-fO2v(indxVoxel,1) fO2v(indxVoxel,1)], ...
        pm.*[-fO2v(indxVoxel,2) fO2v(indxVoxel,2)], ...
        pm.*[-fO2v(indxVoxel,3) fO2v(indxVoxel,3)], ...
        'Color', [0 0 0], 'LineWidth', 2);
end

%% Create DWI

NT.vol = S;
save_nifti(NT, cpDefault.structural_preprocessing.dwiProcessedFile);

NT.vol = NT.vol(:, :, :, 1);
save_nifti(NT, cpDefault.structural_preprocessing.dwiReferenceFile);

%% Create template

% Create segmentationFile
NT.vol = 2*ones(size(S,1),size(S,2),size(S,3)); % White matter
NT.vol(:, :, 1) = 16; % stop region;
NT.vol(:, :, 3) = 16; % stop region;
NT.vol(1:16, 1, 2) = 1001;   NT.vol(1, 1:16, 2) = 1001; % Region 1
NT.vol(17:end, 1, 2) = 1002;   NT.vol(end, 1:16, 2) = 1002; % Region 2
NT.vol(1, 17:end, 2) = 1003;   NT.vol(1:16, end, 2) = 1003; % Region 3
NT.vol(end, 17:end, 2) = 1004;   NT.vol(17:end,end, 2) = 1004; % Region 4
NT.vol(1:3) = 255; % Fake corpus callosum to estimate tensor for CSD

save_nifti(NT, cpDefault.structural_preprocessing.segmentationFile);
save_nifti(NT, strrep(cpDefault.parcellation.parcellationFile, ...
    'TEMPLATE', 'toyAtlas'));

cpDefault.general.templatesDir = 'templateDir';
cpDefault.general.templates = 'toyAtlas';

mkdir(fullfile(cpDefault.general.templatesDir, ...
    cpDefault.general.templates));

ctab = ['1001  region1                          25   5  25    0'; ...
      '1002  region2                          25   5  25    0'; ...
      '1003  region3                          25   5  25    0'; ...
      '1004  region4                          25   5  25    0'];
  
dlmwrite(fullfile(cpDefault.general.templatesDir, ...
    cpDefault.general.templates, 'toyAtlas.annot.ctab'), ctab, '');

dlmwrite(fullfile(cpDefault.general.templatesDir, ...
    cpDefault.general.templates, 'ROIs_toyAtlas.txt'), ...
    [1001; 1002; 1003; 1004], '');

if showFig
x = NT.vol(:,:, 2);
[x2, ~, J] = unique(x);
xx = 1:length(unique(x));
xx = reshape(xx(J), [32 32]);
figure('color', 'white');
hi = imagesc(xx);
cm = cbrewer('qual', 'Accent', 8);
cm = [1 1 1; cm];
colormap(cm);
end

%% Create bvals and bvecs file
dlmwrite(cpDefault.structural_preprocessing.processedBvalsFile, ...
    gtab.bvals);
dlmwrite(cpDefault.structural_preprocessing.processedBvecsFile, ...
    gtab.bvecs);

%% Expected Diffusion peaks
maskCrossing = fO1v(:,1) == fO2v(:, 1);
diffusionPeaks = cat(3, [fO1v; fO1v; fO1v], [fO2v; fO2v; fO2v]);
diffusionPeaks([maskCrossing; maskCrossing; maskCrossing], :, 2) = 0;

for im = 1:length(cpDefault.general.reconstructionMethods)
    m = cpDefault.general.reconstructionMethods{im};
    save(strrep(cpDefault.reconstruction_diffusion.diffusionPeaksFile, ...
        'METHOD', m), 'diffusionPeaks');
end

diffusionPeaks = diffusionPeaks(:, :, 1);
save(strrep(cpDefault.reconstruction_diffusion.diffusionPeaksFile, ...
    'METHOD', 'dti'), 'diffusionPeaks');

%% Create diffusion measures
maskSingle = fO1v(:,1) == fO2v(:, 1) & fO1v(:,1) ~= 0;
maskCrossing = fO1v(:,1) ~= fO1v(:,1);

l1 = -log(SimulateDWMRI([1 0 0],[1 0 0])) / 1500;
l2 = -log(SimulateDWMRI([1 0 0],[0 1 0])) / 1500;
FA = sqrt(1/2) * sqrt((l1 - l2)^2 + (l2 - l1)^2) / sqrt((l1^2 + 2 * l2^2));
RD = l2;
AD = l1;

diffusionMeasures = zeros(size(fO1v,1), 3);
diffusionMeasures(:, 1) = FA * maskSingle + 0.5*FA * maskCrossing;
diffusionMeasures(:, 2) = AD * maskSingle + 0.5*AD * maskCrossing;
diffusionMeasures(:, 3) = RD * maskSingle + 0.5*RD * maskCrossing;

diffusionMeasures = repmat(permute(reshape(diffusionMeasures, ...
    32, 32, []), [1 2 4 3]), [1 1 3 1]);
weightDescriptions = {'fractional anisotropy', ...
            'axial diffusivity', ...
            'radial diffusivity'}';
    
save(cpDefault.reconstruction_diffusion.diffusionMeasuresFile, ...
    'diffusionMeasures', 'weightDescriptions');
