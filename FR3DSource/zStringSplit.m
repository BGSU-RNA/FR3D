% zStringSplit(s,delimiter) splits string s using delimiter, returning a cell array
% Use char(9) for tab

function [parts] = zStringSplit(s,delimiter)

if nargin < 2,
	delimiter = ' ';
end

remain = s;
c = 0;

if isempty(strfind(s,delimiter)),
	parts{1} = s;
else
	while ~isempty(remain),
		c = c + 1;
		[a,remain] = strtok(remain,delimiter);
		parts{c} = a;
	end
end
