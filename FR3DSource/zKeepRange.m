% zKeepZeros(D,E,a,b) returns D with the entries corresponding to entries of E in the range [a,b] remaining

function [D] = zKeepRange(D,E,a,b)

[s,t] = size(D);

[i,j,d] = find(D);                       % non-zero entries of D
[p,q]   = find(E);                       % non-zero entries of E
[c,k]   = setdiff([i j], [p q], 'rows'); % non-zeros of D that are not non-zeros of E
D       = sparse(i(k),j(k),d(k),s,t);    % reconstitute D

