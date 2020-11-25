function S = list2struct(list)
% LIST2STRUCT Convert a paramater name-value cell array to structure.
%
%   S = LIST2STRUCT(LIST) converts a Nx2 cell array to a structure S with 
%   fields named after the parameter names in the first colom of LIST and
%   associated values to the parameter values in the second colom of LIST.
%
%   Input variables:
%   LIST            Nx2 cell array describing N parameter name-value pairs.
%                   The first column is the name of the parameter.
%                   The second column is the value of the parameters.
%
%   Output variables:
%   S               Structure with fields named after the parameter names
%                   and their associated values.

S = struct();
for i = 1:size(list, 1)
    fn = strsplit(list{i, 1}, '.');
    S = setfield(S, fn{:}, list{i, 2}); 
end
