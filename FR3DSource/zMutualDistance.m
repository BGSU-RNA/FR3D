% zMutalDistance(A,L) finds the mutual distances between the rows of A and
% returns a sparse matrix D in which all entries are less than L

function [D] = zMutualDistance(A,L)

if nargin < 2,
  L = Inf;
end

[s,t] = size(A);                       % s is the number of rows

D = zDistance(A);
D = sparse(D .* (D < L));
