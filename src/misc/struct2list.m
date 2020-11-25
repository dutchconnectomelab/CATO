function R = struct2list(S)
% STRUCT2LIST Convert a structure to a paramater name-value cell array.
%
%   R = STRUCT2LIST(S) converts structure S to a Nx2 cell array R with in
%   the first column the field names and in the second column the field
%   values. STRUCT2LIST is recursive and converts also fields of fields.
%   For example:
%       >>test1.test2.test3 = 'hi'
%       >>struct2list(test1)
%             {'test2.test3'}    {'hi'}
%
%   Input variables:
%   S               Structure with N fields.
%
%   Output variables:
%   LIST            Nx2 cell array describing N parameter name-value pairs.
%                   The first column are parameter names. The second column
%                   are parameters values.

R = cell(0,2);
fn1 = fieldnames(S);

for i = 1:length(fn1)
    if isstruct(S.(fn1{i}))
        thisR = struct2list(S.(fn1{i}));
        thisR = [strcat(fn1{i}, '.', thisR(:, 1)) thisR(:, 2)];
        R(end+1:end+(size(thisR, 1)), :) = thisR;
    else
        R(end+1, :) = {fn1{i}, S.(fn1{i})};
    end
    
end
