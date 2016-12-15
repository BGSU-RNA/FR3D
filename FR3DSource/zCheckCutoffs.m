% zCheckCutoffs(D,Normal,Ang,Gap,B) finds the categories whose cutoffs
% include the given x,y displacement D, Normal, angle Ang, and Gap, according to
% the cutoffs in matrix B

function [a] = zCheckCutoffs(D,Normal,Ang,Gap,B)

a = [];                                         % classifications

r = (D(1) > B(:,2)) .* (D(1) < B(:,3)) .* (D(2) > B(:,4)) .* (D(2) < B(:,5));
                                           % check only rows satisfying these
i = find(r);

for j = 1:length(i),
  k = i(j);
  if B(k,10) < B(k,11),
    anglecriterion = (Ang > B(k,10)) & (Ang < B(k,11));
  else
    anglecriterion = (Ang > B(k,10)) || (Ang < B(k,11));
  end

  if ...
    (D(3) > B(k,6)) && ...
    (D(3) < B(k,7)) && ...
    (Normal(3) > B(k,8)) && ...
    (Normal(3) < B(k,9)) && ...
    anglecriterion && ...
    (B(k,12) < 0.1 || abs(Gap) < B(k,12)),    % gap cutoff, if implemented
    a = [a B(k,1)];                          % assign classification in this row of B
  end
end

if length(a) == 0,
  a = 30;                             % no classification given
end
