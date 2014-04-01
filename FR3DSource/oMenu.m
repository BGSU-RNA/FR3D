
function [k] = oMenu(MenuTitle,Buttons,EnterMeans)

if nargin < 3,
	EnterMeans = 1;
end

fprintf('%s\n', MenuTitle);
for i = 1:length(Buttons),
	fprintf('%2d %s\n', i, Buttons{i});
end

zFlushOutput

k = 0;
while k < 1 || k > length(Buttons),
	k = input('Enter the number of your choice here: ','s');
	if isempty(k),
		k = EnterMeans;
	else
		k = str2num(k);
		if isempty(k),
			k = 0;
		end
	end
end
