function [pretty] = prettyjson(ugly)
% PRETTYJSON Makes JSON strings (relatively) pretty.
%
%   INPUT VARIABLES
%   ugly:
%   Char in JSON-format.
%
%   OUTPUT VARIABLES
%   pretty:
%   Prettyfied version of ugly.

%   Inspired by:
%       Yury Bondarenko (2019). prettyjson.m
%       (https://www.github.com/ybnd/prettyjson.m), GitHub. Retrieved
%       February 20, 2020. Shared under the MIT License.

TAB = '   ';

ugly = strrep(ugly, '{', sprintf('{\n'));
ugly = strrep(ugly, '}', sprintf('\n}'));
ugly = strrep(ugly, ',"', sprintf(', \n"'));
ugly = strrep(ugly, ',{', sprintf(', \n{'));
ugly = strrep(ugly, '[', sprintf('[\n'));
ugly = strrep(ugly, ']', sprintf('\n]'));
ugly = strrep(ugly, ',[', sprintf(', \n['));

lines = splitlines(ugly);

balanceItems = 0;


for i = 1:length(lines)
    line = lines{i};
    
    
    % if line starts with closing brac(ket)s
    if ~startsWith(line, {'}', ']'})
    indent = repmat(TAB, [1 balanceItems]);
    end
 
    
    % Prepare for the next line
    % Count brackets and braces
    openItems = length(strfind(line, '[')) + length(strfind(line, '{'));
    closeItems = length(strfind(line, ']')) + length(strfind(line, '}'));
    
    balanceItems = balanceItems + openItems - closeItems;
    
    if startsWith(line, {'}', ']'})
    indent = repmat(TAB, [1 balanceItems]);
    end    
    
    lines{i} = [indent, line];
    
end

pretty = strjoin(lines, newline);
