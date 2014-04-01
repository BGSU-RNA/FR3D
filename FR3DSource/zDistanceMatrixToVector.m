% zSquareform essentially replicates the functionality of Matlab's squareform
% code is taken from http://stackoverflow.com/questions/9659536/alternative-to-using-squareform-matlab

function [v] = zDistanceMatrixToVector(D)

[s,t] = size(D);

if s == t,
	v = zeros(1,(s^2 - s)/2);    % vector to store distances
	c = 1;                       % counter for where to store distances
	for i = 1:(s-1),
		v(c:(c+s-i-1)) = D(i,(i+1):s);
		c = c + s - i;
	end
end
