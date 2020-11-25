function make_structural(catoVersion)
% MAKE_STRUCTURAL   Compile structural pipeline.
%
%   make_structural(CATOVERSION) compiles the structural pipeline for a
%   given version CATOVERSION.

%% Initialization
if ismac
    cleanFiles = {'run_structural_pipeline.sh' 'mccExcludedFiles.log', 'requiredMCRProducts.txt', 'readme.txt'};
    cleanDir = 'structural_pipeline.app';
    osName = 'macos';
    targetFile = 'structural_pipeline.app/Contents/MacOS/structural_pipeline';
    sedFix = '-i ''''';
elseif isunix
    cleanFiles = {'run_structural_pipeline.sh' 'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt'};
    osName = 'unix';
    targetFile = 'structural_pipeline';
    sedFix = '-i';
else
    error('OS is not Unix or Mac');
end

catoDir = fileparts(mfilename('fullpath'));

oldPwd = pwd;
cd(catoDir);
cleanupPwd = onCleanup(@() cd(oldPwd));

targetDir = fullfile(catoDir, 'bin', osName, catoVersion);

[~, ~] = mkdir(targetDir);

versionCompiler = ver;
indx = strcmp({versionCompiler.Name}, 'MATLAB Compiler');
versionCompiler = [versionCompiler(indx).Name, ' ', ...
    versionCompiler(indx).Version, ' ', ...
    versionCompiler(indx).Release];

fid = fopen(fullfile(targetDir, 'VERSION'), 'w+');
fprintf(fid, 'CATO %s\n', catoVersion);
fprintf(fid, 'Operating system: %s (%s)\n', osName, computer);
fprintf(fid, 'Compiled on:      %s\n', datetime('now'));
fprintf(fid, 'Compiled using:   %s', versionCompiler);

copyfile('LICENSE', targetDir);

oldPath = path;
cleanupPath = onCleanup(@() addpath(oldPath));
restoredefaultpath;
addpath(genpath(fullfile(catoDir, 'src')));

%% Compile

mcc('-N', '-d', targetDir, ...
    '-mv', 'src/structural_pipeline/structural_pipeline.m', ...
    '-a', 'src/structural_pipeline', ...
    '-a', 'src/anatomical_steps', ...
    '-a', 'src/misc', ...
    '-a', fullfile(matlabroot, 'toolbox/stats/stats/squareform.m'), ...
    '-R', '-nodisplay', '-R', '-nojvm');
% TODO: Check the number of cores used.

%% Cleanup
movefile(fullfile(targetDir, targetFile), fullfile(targetDir, 'structural_pipeline_bin'));
copyfile('make.template', fullfile(targetDir, 'structural_pipeline'));
system(['sed ' sedFix ' -e ''s/VERSION/', catoVersion, '-', osName, '/g'' '...
    fullfile(targetDir, 'structural_pipeline')]);
system(['cat ' fullfile(targetDir, 'structural_pipeline_bin') ' >> ' fullfile(targetDir, 'structural_pipeline')]);
delete(fullfile(targetDir, 'structural_pipeline_bin'));
system(['chmod +x ' fullfile(targetDir, 'structural_pipeline')]);

cleanFiles = fullfile(targetDir, cleanFiles);
delete(cleanFiles{:})

if ismac
    cleanDir = fullfile(targetDir, cleanDir);
    rmdir(cleanDir, 's');
elseif isunix
    delete(fullfile(targetDir, 'functional_pipeline_bin'));
end

mkdir(fullfile(targetDir, 'templates'));
copyfile(fullfile(catoDir, 'src', 'templates'), fullfile(targetDir, 'templates'));
mkdir(fullfile(targetDir, 'structural_preprocessing'));
copyfile(fullfile(catoDir, 'src', 'structural_preprocessing'), fullfile(targetDir, 'structural_preprocessing'));

zip(fullfile(catoDir, 'bin', osName, ['CATO-' catoVersion '-' osName '.zip']), fullfile(targetDir, '*'));

fprintf('done.\n');