function updatedParams = updateParameters(oldParams, newParams)
% UPDATEPARAMETERS Update OLDPARAMS with NEWPARAMS.
%
%   updatedParams = updateParameters(oldParams, newParams) adds newParams
%   parameters to oldParams and saves them in updatedParams. Parameters
%   that already exists in oldParams are overwritten by newParams.
%
%   INPUT VARIABLES
%   oldParams:
%   Old parameters that are updated.
%
%   newParams:
%   New parameters that are added to (or overwrite) parameters in
%   oldParams.
%
%   OUTPUT VARIABLES
%   updatedParams:
%   Struct with upddated parameters.

if ~isstruct(oldParams)
    oldParams = list2struct(oldParams);
end

if isstruct(newParams)
    newParams = struct2list(newParams);
end

updatedParams = oldParams;

for iParam = 1:size(newParams, 1)
    fn = strsplit(newParams{iParam, 1}, '.');
    updatedParams = setfield(updatedParams, fn{:}, newParams{iParam, 2});
end

end