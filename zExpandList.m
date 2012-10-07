% zExpandList(Indices,n) lists all indices within n of Indices

function [NewIndices] = zExpandList(Indices,n,IMax)

s = sort(Indices);

if s(1) > 1,
  NewIndices = max(1,(s(1)-n)):(s(1)-1);
else
  NewIndices = [];
end

for i = 2:length(s),
  if s(i)-s(i-1) > 1,
    NewIndices = [NewIndices (s(i-1)+1):min(s(i-1)+n,s(i)-n-1) max(s(i)-n,s(i-1)+1):(s(i)-1)];
  end
end

if s(end) < IMax,
  NewIndices = [NewIndices (s(end)+1):min(IMax,s(end)+n)];
end
