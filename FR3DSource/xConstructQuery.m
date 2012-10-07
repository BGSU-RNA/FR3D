% xConstructQuery(Query,File) fills in details of Query from File

% defaults:
% Query.Geometric = 1 if Query.Filename and Query.NTList are defined
% Query.Geometric = 0 if not; in that case, Query.Mask, Query.Edges, or
%      Query.MaxDiff need to be defined, or there is no focus for the search
% Query.ExcludeOverlap = 1 if there are 7 or more nucleotides
% Query.Mask is "NNNNN..."
% Query.Sequential is 1 if Query.MaxDiff is non-empty
% Query.MaxDiff is [Inf Inf ...] by default
% Query.LocWeight is [1 1 1 1...]
% Query.AngleWeight is [1 1 1 1...]
% Query.DiscCutoff is 0.4
% Query.RelCutoff is Query.DiscCutoff, if not overridden
% Query.ChainList is not needed unless there is ambiguity in nucleotide
%                 numbers, and if there is, xConstructQuery will tell you

function [Query] = xConstructQuery(Query,File)

% ------------------------------- Determine whether search is geom. or symb.

if ~isfield(Query,'Geometric'),
  if isfield(Query,'Filename') & isfield(Query,'NTList'),
    Query.Geometric = 1;
  else
    Query.Geometric = 0;
  end
end

% ------------------------------- Determine number of nucleotides if geometric

if Query.Geometric == 1,
  Query.NumNT = length(Query.NTList);
end  

% ------------------------------- Determine number of nucleotides if symbolic

if Query.Geometric == 0,
  m = 0;
  if isfield(Query,'Edges'),
    [s,t] = size(Query.Edges);
    m = max([m s t]);               % number of nucleotides
  end
  if isfield(Query,'Diagonal'),
    m = length(Query.Diagonal);
  end
  if isfield(Query,'Mask'),
    m = length(Query.Mask);
  end
  if isfield(Query,'Diff'),
    [s,t] = size(Query.Diff);
    m = max([m s t]);               % number of nucleotides
  end
  if isfield(Query,'MaxDiff'),
    [s,t] = size(Query.MaxDiff);
    m = max([m s t]);
  end
  if isfield(Query,'MaxDiffVector'),
    [s,t] = size(Query.MaxDiffVector);
    m = max(m,max(s,t)+1);
  end
  if isfield(Query,'MinDiff'),
    [s,t] = size(Query.MinDiff);
    m = max([m s t]);
  end
  if isfield(Query,'MinDiffVector'),
    [s,t] = size(Query.MinDiffVector);
    m = max(m,max(s,t)+1);
  end
  if m == 0,
    fprintf('Not enough information to form a query motif.\n');
    fprintf('No search was conducted.\n');
  else
    Query.NumNT = m;
  end
end

if isfield(Query,'NumNT'),

% --------- Set default values

if ~isfield(Query,'ExcludeOverlap'),
  if Query.NumNT >= 7,
    Query.ExcludeOverlap = 1;
  else
    Query.ExcludeOverlap = 0;
  end
end

if ~isfield(Query,'Number'),
  Query.Number = '0';
end

if ~isfield(Query,'Name'),
  Query.Name = num2str(Query.Number);
end

if ~isfield(Query,'Description'),
  Query.Description = Query.Name;
end

if ~isfield(Query,'Config'),
  for i=1:Query.NumNT,
    Query.Config{i} = '';
  end
else
  for i=(length(Query.Config)+1):Query.NumNT,
    Query.Config{i} = '';
  end
end

if (Query.Geometric == 0) & ~isfield(Query,'Diameter'),
  Query.Diameter = 30;
end

% ---------------------------------------------------------------------

if Query.Geometric == 1,                    % model comes from a file
    
  if ~isfield(Query,'Diagonal'),
    if ~isfield(Query,'LocWeight'),
      Query.LocWeight = ones(1,Query.NumNT);
    end
    if ~isfield(Query,'AngleWeight'),
      Query.AngleWeight = ones(1,Query.NumNT);
    end
  end

  if ~isfield(Query,'DiscCutoff'),
    Query.DiscCutoff = 0.4;
  end

  if ~isfield(Query,'RelCutoff'),
    Query.RelCutoff = Query.DiscCutoff;
  end

  % --------- basic parameters for the model

  if isfield(Query,'ChainList'),
    Query.Indices  = zIndexLookup(File,Query.NTList,Query.ChainList);
  else
    Query.Indices  = zIndexLookup(File,Query.NTList);
  end

  for i=1:Query.NumNT,
    Query.NT(i) = File.NT(Query.Indices(i));
  end

  Query.Edge  = File.Edge(Query.Indices,Query.Indices);
end

% --------- Interpret model mask

% For nucleotide i, the vector Query.OKCodes{i} has 0's and 1's 
% which tell which nucleotide codes (A=1, C=2,...) meet the mask

if isfield(Query,'Mask'),
  for i=1:Query.NumNT,
    switch Query.Mask(i)
     case 'A', OK = [1 0 0 0];
     case 'C', OK = [0 1 0 0];
     case 'G', OK = [0 0 1 0];
     case 'U', OK = [0 0 0 1];
     case 'M', OK = [1 1 0 0];
     case 'R', OK = [1 0 1 0];
     case 'W', OK = [1 0 0 1];
     case 'S', OK = [0 1 1 0];
     case 'Y', OK = [0 1 0 1];
     case 'K', OK = [0 0 1 1];
     case 'V', OK = [1 1 1 0];
     case 'H', OK = [1 1 0 1];
     case 'D', OK = [1 0 1 1];
     case 'B', OK = [0 1 1 1];
     case 'N', OK = [1 1 1 1];
    end
    Query.OKCodes{i} = OK;
  end
end

% --------- Interpret model diagonal (mask and weights)

% For nucleotide i, the vector Query.OKCodes{i} has 0's and 1's 
% which tell which nucleotide codes (A=1, C=2,...) meet the mask

if isfield(Query,'Diagonal'),
  for i=1:Query.NumNT,
    [OK,LW,AW]           = xGetNuclSpec(Query.Diagonal{i});
    Query.OKCodes{i}     = OK;
    Query.LocWeight(i)   = LW;
    Query.AngleWeight(i) = AW;
  end
end

% --------- If no nucleotide mask has been encountered

if ~isfield(Query,'OKCodes'),
  for i=1:Query.NumNT,
    Query.OKCodes{i} = [1 1 1 1];
  end
end

% --------- Make edge matrix antisymmetric and numeric

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

RevPair = [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];  % exchange AC for CA etc.

if isfield(Query,'Edges'),
  Query.Edges{Query.NumNT,Query.NumNT} = 0;
  Query.EdgeNums{Query.NumNT,Query.NumNT} = 0;
  Query.ExcludeEdges{Query.NumNT,Query.NumNT} = 0;
  Query.OKPairs{Query.NumNT,Query.NumNT} = 0;
  Query.ExPairs{Query.NumNT,Query.NumNT} = 0;
  Query.BasePhos{Query.NumNT,Query.NumNT} = 0;
  Query.ExcludeBasePhos{Query.NumNT,Query.NumNT} = 0;
  Query.Flank{Query.NumNT,Query.NumNT} = [];
  Query.Range{Query.NumNT,Query.NumNT} = [];
  Query.Coplanar{Query.NumNT,Query.NumNT} = [];
  Query.OKBB{Query.NumNT,Query.NumNT} = [];
  Query.ExBB{Query.NumNT,Query.NumNT} = [];

  for i=1:Query.NumNT,
    for j=(i+1):Query.NumNT,
      if ~isempty(Query.Edges{j,i}),
        [ReqEdge,ExEdge,OKPairs,ExPairs,BP1,BP2,EBP1,EBP2,Flank,Range,Coplanar,ReqBB,ExBB] = xGetEdgeNums(Query.Edges{j,i});
        Query.EdgeNums{j,i}     =  ReqEdge;
        Query.EdgeNums{i,j}     = -ReqEdge;
        Query.ExcludeEdges{j,i} =  ExEdge;
        Query.ExcludeEdges{i,j} = -ExEdge;
        Query.OKPairs{j,i}      = OKPairs;
        Query.OKPairs{i,j}      = RevPair(OKPairs);
        Query.ExPairs{j,i}      = ExPairs;
        Query.ExPairs{i,j}      = RevPair(ExPairs);
        Query.BasePhos{j,i}     = BP1;
        Query.BasePhos{i,j}     = BP2;
        Query.ExcludeBasePhos{j,i} = EBP1;
        Query.ExcludeBasePhos{i,j} = EBP2;
	Query.Flank{i,j}           = Flank;
	Query.Range{i,j}           = Range;
        Query.Coplanar{i,j}        = Coplanar;
        Query.OKBB{j,i}      = ReqBB;
        Query.ExBB{j,i}      = ExBB;
        Query.OKBB{i,j}      = ReqBB;
        Query.ExBB{i,j}      = ExBB;
      elseif ~isempty(Query.Edges{i,j}),
        [ReqEdge,ExEdge,OKPairs,ExPairs,BP1,BP2,EBP1,EBP2,Flank,Range,Coplanar,ReqBB,ExBB] = xGetEdgeNums(Query.Edges{i,j});
        Query.EdgeNums{j,i}     = -ReqEdge;
        Query.EdgeNums{i,j}     =  ReqEdge;
        Query.ExcludeEdges{j,i} = -ExEdge;
        Query.ExcludeEdges{i,j} =  ExEdge;
        Query.OKPairs{j,i}      = RevPair(OKPairs);
        Query.OKPairs{i,j}      = OKPairs;
        Query.ExPairs{j,i}      = RevPair(ExPairs);
        Query.ExPairs{i,j}      = ExPairs;
        Query.BasePhos{i,j}     = BP1;
        Query.BasePhos{j,i}     = BP2;
        Query.ExcludeBasePhos{i,j} = EBP1;
        Query.ExcludeBasePhos{j,i} = EBP2;
        Query.Flank{j,i}           = Flank;
        Query.Range{j,i}           = Range;
        Query.Coplanar{j,i}        = Coplanar;
        Query.OKBB{i,j}      = ReqBB;
        Query.ExBB{i,j}      = ExBB;
        Query.OKBB{j,i}      = ReqBB;
        Query.ExBB{j,i}      = ExBB;
      else
        Query.EdgeNums{i,j} = [];
        Query.EdgeNums{i,j} = [];
      end
    end
  end
end

% --------- precompute parameters for geometric screening and ranking

if Query.Geometric > 0,
  Query.LocWeight = Query.NumNT * Query.LocWeight / sum(Query.LocWeight);
                                        % weights sum to Query.NumNT
  Query.SSCutoff  =(Query.NumNT^2)*(Query.RelCutoff^2)*cumsum(Query.LocWeight);
                                  % cutoffs for first 1, 2, 3, ... nucleotides
  
  Query = xPrecomputeForDiscrepancy(Query);

  Query.Distance = zMutualDistance(Query.Centers,Inf);

  Query.LDiscCutoff = (Query.NumNT*Query.RelCutoff)^2;

  if (Query.NumNT >= 2),

    if isfield(Query,'Flex'),                  % not checked in a long time
      Query.DistanceScreen = 0;        % cannot use sums of squares with flex
      Query.Flex(Query.NumNT,Query.NumNT) = 0; % make big enough
      Query.Flex = Query.Flex + Query.Flex';   % make symmetric
    end

    Query.DistCutoff = max(max(Query.Distance)) ...
                     + sqrt(2)*Query.NumNT*Query.DiscCutoff;
                                     % distances this large needed in File

  else       % --------- Special calculations for two-nucleotide motifs
    Query.DistCutoff = Query.Distance(1,2) + 2 * Query.DiscCutoff;
  end

else
  Query.SSCutoff = Inf * ones(1,Query.NumNT);
  Query.DistCutoff = 30;        % default distance cutoff for symbolic search
end

% --------- Read minimum and maximum distance specifications

if isfield(Query,'Diff'),
  [s,t] = size(Query.Diff);
  if ~strcmp(class(Query.Diff),'char'),
    for i=1:s,
      for j=1:t,
        if ~isempty(Query.Diff{i,j}),
          [md,MD,sign] = xGetDiffSpec(Query.Diff{i,j});
          Query.MaxDiff(i,j) = MD;
          Query.MinDiff(i,j) = md;
          Query.DifferenceSign(i,j) = sign;
        end
      end
    end
  end
end

% --------- implement maxdiff screening and set difference signs

if isfield(Query,'MaxDiff'),               % matrix of max differences
  D  = Inf*ones(Query.NumNT,Query.NumNT);
  DS = zeros(Query.NumNT,Query.NumNT);

  [s,t] = size(Query.MaxDiff);

    for i=1:s,
      for j=1:t,
        D(i,j) = Query.MaxDiff(i,j);       % read in values set
        if D(i,j) == 0,                    % if no value set, use infinity
          D(i,j) = Inf;
        end
        DS(i,j) = Query.DifferenceSign(i,j);
     end
    end

    for i=1:Query.NumNT,
      D(i,i) = 0;
      for j=1:Query.NumNT,
        D(i,j) = min(D(i,j),D(j,i));       % symmetrize entries in D
        D(j,i) = D(i,j);

        if DS(i,j) ~= 0,
          DS(j,i) = -DS(i,j);
        end
      end
    end

   for m=1:Query.NumNT,                    % just in case
    for i=1:Query.NumNT,
      for j=1:Query.NumNT,
        for k=1:Query.NumNT,
          D(i,j) = min(D(i,j),D(i,k)+D(k,j)); % infer shorter maxdiff's
        end
      end
    end
   end

   Query.MaxDiffMat = D;
   Query.DifferenceSignMat = DS;
end

% --------- implement mindiff screening

if isfield(Query,'MinDiff'),
  D = zeros(Query.NumNT,Query.NumNT);       % minimum distance 1

  [s,t] = size(Query.MinDiff);

    for i=1:s,
      for j=1:t,
        D(i,j) = Query.MinDiff(i,j);       % read in values set
        if D(i,j) == 0,                    % if no value set, use 1
          D(i,j) = 1;
        end
     end
    end

    for i=1:Query.NumNT,
      D(i,i) = 0;
      for j=1:Query.NumNT,
        D(i,j) = max(D(i,j),D(j,i));       % symmetrize entries in D
        D(j,i) = D(i,j);
      end
    end

    Query.MinDiffMat = D;
end

% --------- implement maxdiff screening for a maxdiff vector

if isfield(Query,'MaxDiffVector'),
  D = Inf*ones(Query.NumNT,Query.NumNT);

  [s,t] = size(Query.MaxDiffVector);

  if min(s,t) == 1,                        % vector of max differences

    for i=1:(Query.NumNT-1),
      D(i,i+1) = Query.MaxDiffVector(i);   % put values in a matrix
    end

    for d=3:Query.NumNT,                   % diagonal to consider
      for i=1:(Query.NumNT-d+1),
        j = i + d - 1;
        D(i,j) = D(i,j-1) + D(i+1,j);      % infer max differences
      end
    end

    Query.MaxDiffMat = triu(D) + triu(D)';    % make symmetric

    if isfield(Query,'Filename'),
      [y,i] = sort(Query.Indices);         % re-order for nucleotide order
    else
      i = 1:Query.NumNT;
    end

    k(i) = 1:Query.NumNT;

    Query.MaxDiffMat = Query.MaxDiffMat(k,k); % re-order matrix of differences
  end    
end

% --------- implement mindiff screening for a mindiff vector

if isfield(Query,'MinDiffVector'),
  D = zeros(Query.NumNT,Query.NumNT);       % minimum distance 1

  [s,t] = size(Query.MinDiffVector);

  if min(s,t) == 1,                         % vector of min differences

    for i=1:(Query.NumNT-1),
      D(i,i+1) = Query.MinDiffVector(i);    % put values in a matrix
    end

    Query.MinDiffMat = D + D';              % make symmetric

    if isfield(Query,'Filename'),
      [y,i] = sort(Query.Indices);          % re-order for nucleotide order
    else
      i = 1:Query.NumNT;
    end

    k(i) = 1:Query.NumNT;

    Query.MinDiffMat = Query.MinDiffMat(k,k); % re-order matrix of differences
  end    
end

end  % if isfield(Query,NumNT)
