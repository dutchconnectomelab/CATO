function updateStatus(statusFile, status)
% UPDATESTATUS Update status file with current pipeline status.
%
%   INPUT VARIABLES
%   statusFile:
%   Status file with extension .running, .error or .finished (in JSON-format).
%
%   status:
%   Struct with current pipeline status.

statusFileOld = {strrep(statusFile, 'STATUS', 'running'), ...
    strrep(statusFile, 'STATUS', 'error'), ...
    strrep(statusFile, 'STATUS', 'finished')};
statusFileOld = statusFileOld(isfile(statusFileOld));

if length(statusFileOld) == 1
    % status file already exists
    statusUpdated = readConfigFile(statusFileOld{1});
    statusUpdated = updateParameters(statusUpdated, status);
    
    delete(statusFileOld{1})
    
elseif isempty(statusFileOld)
    statusUpdated = status;
else
    error('CATO:updateStatus:MultipleLogFiles', ...
        'Multiple log files found. \n %s', strjoin(statusFileOld));
end

statusFileNew = strrep(statusFile, 'STATUS', statusUpdated.general);

try
    saveConfigFile(statusFileNew, statusUpdated)
catch
    error('CATO:updateStatus:CouldNotOpenLogFile', ...
        'Could not create status file (%s)', statusFileNew)
end





