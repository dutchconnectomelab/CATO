function make_functional(catoVersion)
% MAKE_FUNCTIONAL   Compile functional pipeline.
%
%   make_functional(CATOVERSION) compiles the functional pipeline for a
%   given version CATOVERSION.

%% Initialization
if ismac
    cleanFiles = {'run_functional_pipeline.sh' 'mccExcludedFiles.log', 'requiredMCRProducts.txt', 'readme.txt'};
    cleanDir = 'functional_pipeline.app';
    osName = 'macos';
    targetFile = 'functional_pipeline.app/Contents/MacOS/functional_pipeline';
    sedFix = '-i ''''';
elseif isunix
    cleanFiles = {'run_functional_pipeline.sh' 'mccExcludedFiles.log' 'readme.txt' 'requiredMCRProducts.txt'};    
    osName = 'unix';
    targetFile = 'functional_pipeline';
    sedFix = '-i';
else
    error('OS is not unix or Mac');
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
restoredefaultpath;
cleanupPath = onCleanup(@() addpath(oldPath));


%% Compile

mcc('-N', '-d', targetDir, ...
    '-mv', 'src/functional_pipeline/functional_pipeline.m', ...
    '-a', 'src/functional_pipeline', ...    
    '-a', fullfile(matlabroot, 'toolbox/stats/stats/squareform.m'), ...
    '-a', fullfile(matlabroot, 'toolbox/signal/signal/butter.m'), ...
    '-a', fullfile(matlabroot, 'toolbox/signal/signal/filtfilt.m'), ...
    '-a', 'src/anatomical_steps', ...
    '-a', 'src/misc', ...    
    '-R', '-nodisplay', '-R', '-nojvm');
% TODO: Check the number of cores used.

%% Cleanup

movefile(fullfile(targetDir, targetFile), fullfile(targetDir, 'functional_pipeline_bin'));
copyfile('make.template', fullfile(targetDir, 'functional_pipeline'));
system(['sed ' sedFix ' -e ''s/VERSION/', catoVersion, '-', osName, '/g'' '...
    fullfile(targetDir, 'functional_pipeline')]);
system(['sed ' sedFix ' -e ''s/structural/functional/g'' '...
    fullfile(targetDir, 'functional_pipeline')]);
system(['cat ' fullfile(targetDir, 'functional_pipeline_bin') ' >> ' fullfile(targetDir, 'functional_pipeline')]);
delete(fullfile(targetDir, 'functional_pipeline_bin'));
system(['chmod +x ' fullfile(targetDir, 'functional_pipeline')]);

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
mkdir(fullfile(targetDir, 'functional_preprocessing'));
copyfile(fullfile(catoDir, 'src', 'functional_preprocessing'), fullfile(targetDir, 'functional_preprocessing'));

zip(fullfile(catoDir, 'bin', osName, ['CATO-' catoVersion '-' osName '.zip']), fullfile(targetDir, '*'));

fprintf('done.\n');