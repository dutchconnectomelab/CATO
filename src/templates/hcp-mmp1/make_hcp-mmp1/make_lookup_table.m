templateDir = fileparts(fileparts(mfilename('fullpath')));

%% Adjust RAW hcp-mmp1.annot file
% Rename the regions such that they are in line with other atlases and
% such that the function "matchROIs" works in the parcellation step.
%
% RAW data is obtained from: 
% Mills, Kathryn (2016): HCP-MMP1.0 projected on fsaverage. figshare.
% Dataset. https://doi.org/10.6084/m9.figshare.3498446.v2

for hemi = {'lh', 'rh'}

[V, L, ct] = read_annotation(['RAW/' hemi{1} '.HCP-MMP1.annot']);

names = ct.struct_names;
names = strrep(names, '???', 'unkown');
names = regexprep(names, '^L_', '');
names = regexprep(names, '^R_', '');
names = erase(names, '_ROI');
ct.struct_names = names;

write_annotation(fullfile(templateDir, [hemi{1} '.hcp-mmp1.annot']), V, L, ct);

% Check vertices and labels are not changed
[V2, L2, ct2] = read_annotation(fullfile(templateDir, [hemi{1} '.hcp-mmp1.annot']));
assert(isequal(V, V2));
assert(isequal(L, L2));

fprintf('New annot file made for %s.\n', hemi{1})

end

%% Make annot.ctab file
% Merge LUT table from lh and rh annot file and FreeSurferColorLUT
lutFSFile = '/Applications/freesurfer/6.0.0/FreeSurferColorLUT.txt';
lutFile = fullfile(templateDir, 'hcp-mmp1.annot.ctab'); 

[~, ~, cL] = read_annotation(fullfile(templateDir, 'lh.hcp-mmp1.annot'));

namesL = cL.struct_names;
namesL = strcat('ctx-lh-', namesL);
rgbvL = cL.table(:, 1:4);
codesL = (1:length(namesL))'-1;

[~, ~, cR] = read_annotation(fullfile(templateDir, 'rh.hcp-mmp1.annot'));
namesR = cR.struct_names;
namesR = strcat('ctx-rh-', namesR);
rgbvR = cR.table(:, 1:4);
codesR = (1:length(namesR))'-1;

[codesFS, namesFS, rgbvFS] = read_fscolorlut(lutFSFile);
namesFS = mat2cell(namesFS, ones(size(namesFS,1),1), size(namesFS,2));
indxFS = codesFS < 1000;

% Comnbine FS standard (subcortical etc.) with cortical parcellation
codes = [codesFS(indxFS); 1000 + codesL; 2000 + codesR];
names = [namesFS(indxFS); namesL; namesR];
rgbv = [rgbvFS(indxFS, :); rgbvL; rgbvR];

% Write color table
names = pad(names, max(cellfun(@(x) length(x), names)));
fid = fopen(lutFile, 'w');
for i = 1:length(codes)
    fprintf(fid, '%u\t%s\t%u\t%u\t%u\t%u\n', ...
        codes(i), names{i}, rgbv(i, :));
end
fclose(fid);

