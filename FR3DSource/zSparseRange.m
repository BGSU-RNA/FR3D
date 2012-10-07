% zSparseRange returns a sparse matrix A with the values of M satisfying a <= m <= b.  No problem if a = 0. If v = 1, it returns 1's instead of M with values.

function [A] = zSparseRange(M,a,b,v)

if nargin < 4,
  v = 0;
end

[i,j,s] = find(M);
[x,y] = size(M);
k = find( (s >= a) .* (s <= b) );

if v == 0,
  A = sparse(i(k),j(k),s(k),x,y);       % return values of M
else
  A = sparse(i(k),j(k),1,x,y);          % return 1's where non-zero
end
