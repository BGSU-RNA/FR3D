% zSparseValues returns a sparse matrix A with the values m of M satisfying any(m == c)

function [A] = zSparseRange(M,c,v)

if nargin < 3,
  v = 0;
end

[i,j,s] = find(M);
[x,y] = size(M);

OK = zeros(size(s));
for cc = 1:length(c),
  OK = OK + (s == c(cc));
end

k = find( OK );

if v == 0,
  A = sparse(i(k),j(k),s(k),x,y);       % return values of M
else
  A = sparse(i(k),j(k),1,x,y);          % return 1's where non-zero
end
