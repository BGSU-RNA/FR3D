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
elseif delimiter ~= char(9),
	while ~isempty(remain),
		c = c + 1;
		[a,remain] = strtok(remain,delimiter);
		parts{c} = a;
	end
else
	j = find(s == char(9));
	parts{1} = s(1:(j(1)-1));
	for k = 1:(length(j)-1),
		parts{k+1} = s((j(k)+1):(j(k+1)-1));
	end
	parts{length(j)+1} = s((j(end)+1):end);
end
