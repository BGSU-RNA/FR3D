% zOrderbySimilarity uses the result of a cluster analysis and puts the
% instances in order, consistent with the tree, so that similar branches are
% near each other in the list

% DD is a symmetric mutual distance matrix
% method is optional, it is the hierarchical clustering method to use
% W is optional, it is the weight matrix to use
% Verbose is optional, it tells whether or not to display output
% the result minperm is a permutation of the numbers 1:length(DD(:,1))

function [minperm,minscore,mingroup] = zOrderbySimilarity(DD,method,W,Verbose)

[s,t] = size(DD);

if nargin < 4,
  Verbose = 0;
end

if nargin < 3,
  W = 2.^(-abs(ones(s,1)*(1:s) - (1:s)'*ones(1,s)));  % weight matrix
elseif isempty(W),
  W = 2.^(-abs(ones(s,1)*(1:s) - (1:s)'*ones(1,s)));  % weight matrix
end

if nargin < 2,
  method = {'single','average','complete','Ward','weighted','centroid','median'};
elseif isempty(method),
  method = {'single','average','complete','Ward','weighted','centroid','median'};
end

if strcmp(class(method),'char'),
  method = {method};
end

minscore  = Inf;
minmethod = 1;
minperm   = 1:s;
mingroup  = {};

if s > 1,

for m = 1:length(method),

  % ----------------------------------------- Cluster analysis

  D = full(DD);                                  % not sparse, store all zeros
  for i = 1:s,
    D(i,i) = 0;                                  % set diagonal to zero
  end

  Y = squareform(full(D));                       % convert to a vector

  warning off
  Z = linkage(Y,method{m});                      % compute cluster tree

  warning on

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
    g  = Group{Z(z,1)};                   % members of one group
    gg = fliplr(g);                       % reverse order
    h  = Group{Z(z,2)};                   % members of the group being merged

    if ~isempty(intersect(g,h))           % this should never happen!!!
      if Verbose > 0,
        fprintf('Clustering method %10s has non-disjoint groups!\n', method{m});
      end
      h = setdiff(h,g);                   % take out the intersection
    end

    S(1) = Score(D([g h],[g h]),W);       % score four different orderings
    S(2) = Score(D([h g],[h g]),W);
    S(3) = Score(D([gg h],[gg h]),W);
    S(4) = Score(D([h gg],[h gg]),W);

    [y,i] = sort(S);                      % find the lowest score

    switch i(1)                           % merge according to lowest score
      case 1, Group{gc+1} = [g h];
      case 2, Group{gc+1} = [h g];
      case 3, Group{gc+1} = [gg h];
      case 4, Group{gc+1} = [h gg];
    end

    gc = gc + 1;                          % move on to the next group
  end

  q = Group{end};

  S = Score(D(q,q),W);

  % ---------- I do not understand why, but sometimes this procedure
  % ---------- gives a permutation which is too short.

  if length(q) < s,
    if Verbose > 0,
      fprintf('Clustering method %10s misses some rows!\n', method{m});
    end
    S = Inf;                                  % avoid this permutation!
  end
    
  if Verbose > 0,
    fprintf('Clustering method %10s gives score %8.4f\n', method{m}, S);
  end

  if S < minscore,
    minscore  = S;
    minmethod = m;
    minperm   = q;
    mingroup  = Group;
  end    

end

end

if Verbose > 0,
  fprintf('Using method %10s.\n', method{minmethod});
end
  
% -----------------------------------------------------------------------
% The Score function tells which type of distance matrix is preferred.
% This simple scoring scheme prefers small values near the diagonal.

function [S] = Score(M,W)

[s,t] = size(M);

S = sum(sum(M .* W(1:s,1:t)));


