% zOrderbySimilarity uses the result of a cluster analysis and puts the
% instances in order, consistent with the tree, so that similar branches are
% near each other in the list

% DD is a symmetric mutual distance matrix
% q is a permutation

function [q] = zOrderbySimilarity(DD)

[s,t] = size(DD);

W = 2.^(-abs(ones(s,1)*(1:s) - (1:s)'*ones(1,s)));  % weight matrix

% ----------------------------------------- Cluster analysis

D = full(DD);
for i = 1:s,
  D(i,i) = 0;                                  % set diagonal to zero
end

Y = squareform(full(D));                       % convert to a vector
Z = linkage(Y,'average');                      % compute cluster tree

% D is the square matrix giving mutual distances, with zero down the diagonal
% Y is a vector, needed for the linkage function
% Z is the linkage, telling what coalesces with what

% ----------------------------------------- Re-order the elements in groups

q = 1:s;                                % simplest ordering, doesn't do well

for i = 1:s,
  Group{i} = i;                         % initial groups have one element
end

gc = s;                                 % group counter

for z = 1:length(Z(:,1)),
  g  = Group{Z(z,1)};
  h  = Group{Z(z,2)};
  gg = fliplr(g);

  S(1) = Score(D([g h],[g h]),W);
  S(2) = Score(D([h g],[h g]),W);
  S(3) = Score(D([gg h],[gg h]),W);
  S(4) = Score(D([h gg],[h gg]),W);

  [y,i] = sort(S);

  switch i(1)
    case 1, Group{gc+1} = [g h];
    case 2, Group{gc+1} = [h g];
    case 3, Group{gc+1} = [gg h];
    case 4, Group{gc+1} = [h gg];
  end

  gc = gc + 1;
end

q = Group{end};

% -----------------------------------------------------------------------
% The Score function tells which type of distance matrix is preferred.
% This simple scoring scheme prefers small values near the diagonal.

function [S] = Score(M,W)

% S = 2*sum(diag(M,1)) + sum(diag(M,2));

[s,t] = size(M);
S = sum(sum(M .* W(1:s,1:t)));

