% xPairwiseScreen returns a sparse matrix with non-zero entries corresponding to pairs of bases which satisfy all given constraints
% Codes are the base codes from File
% 

function [Screen] = xPairwiseScreen(File,Codes,Query,p,q,PC);

if Query.Geometric == 0,
  D = File.Distance .* (File.Distance < Query.Diameter);
                                        % cap distance for non-geometric search
else
  D = File.Distance;
end

% --------- Screen according to interaction between nucleotides

if isfield(Query,'EdgeNums'),           % if screening by edges, incl stacking
  if length(Query.EdgeNums{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.EdgeNums{p,q}),
      E = E + (fix(File.Edge) == Query.EdgeNums{p,q}(i));
    end
    D = D .* (E > 0);                         % include only those that match
  end
end
    
if isfield(Query,'ExcludeEdges'),                 % if excluding by edges
  if length(Query.ExcludeEdges{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeEdges{p,q}),
      E = E + (fix(File.Edge) == Query.ExcludeEdges{p,q}(i));
    end
    D = D .* (E == 0);
  end
end
    
if isfield(Query,'OKPairs'),                 % if screening by paircode
  if length(Query.OKPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.OKPairs{p,q}),
      E = E + sparse(PC == Query.OKPairs{p,q}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Query,'ExPairs'),                 % if screening by paircode
  if length(Query.ExPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExPairs{p,q}),
      E = E + sparse(PC == Query.ExPairs{p,q}(i));
    end
    D = D .* (E == 0);
  end
end

if isfield(Query,'OKBB'),           % if screening by backbone
  if length(Query.OKBB{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.OKBB{p,q}),
      E = E + (fix(File.Backbone) + fix(File.Backbone') == Query.OKBB{p,q}(i));
    end
    D = D .* (E > 0);                         % include only those that match
  end
end
    
if isfield(Query,'ExBB'),                 % if excluding by backbone
  if length(Query.ExBB{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExBB{p,q}),
      E = E + (fix(File.Backbone) + fix(File.Backbone') == Query.ExBB{p,q}(i));
    end
    D = D .* (E == 0);
  end
end

if isfield(Query,'BasePhos'),                 % if screening by base-phosphate
  if length(Query.BasePhos{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.BasePhos{p,q}),
      E = E + (fix(File.BasePhosphate) == Query.BasePhos{p,q}(i));
    end
    D = D .* (E > 0);
  end
  if length(Query.BasePhos{q,p} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.BasePhos{q,p}),
      E = E + (fix(File.BasePhosphate') == Query.BasePhos{q,p}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Query,'ExcludeBasePhos'),                 % if excluding by edges
  if length(Query.ExcludeBasePhos{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeBasePhos{p,q}),
      E = E + (fix(File.BasePhosphate) == Query.ExcludeBasePhos{p,q}(i));
    end
    D = D .* (E == 0);
  end
  if length(Query.ExcludeBasePhos{q,p} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Query.ExcludeBasePhos{q,p}),
      E = E + (fix(File.BasePhosphate') == Query.ExcludeBasePhos{q,p}(i));
    end
    D = D .* (E == 0);
  end
end

if isfield(Query,'Flank'),
  if ~isempty(Query.Flank{p,q}),
    D = D .* File.Flank;                 % only keep pairs between nested cWW
  end
end

if isfield(Query,'Range'),               % screen by range restriction
  if ~isempty(Query.Range{p,q}),
    R = Query.Range{p,q};
    D = D .* (File.Range >= R(1)) .* (File.Range <= R(2));
  end
end

if isfield(Query,'Coplanar'),            % screen by coplanarity
 if ~isempty(Query.Coplanar{p,q}),

  switch Query.Coplanar{p,q}
  case 0,
    D = D .* (File.Coplanar < 0.5);      % ~cp
  case 1,
    D = D .* (File.Coplanar >= 0.5);     % cp
  case 2,
    D = D .* (File.Coplanar > 0);        % ncp cp
  case 3,
    D = D .* (File.Coplanar > 0) .* (File.Coplanar < 0.5);  % ncp
  case 4,
    D = D .* ((File.Coplanar == 0) + (File.Coplanar >= 0.5)); % ~ncp
  case 5,
    D = D .* (File.Coplanar == 0);       % ~ncp ~cp
  end

%figure(3)
%hist(nonzeros(D),30)

 end
end

    
[i,j] = find(D);         % nucleotide pairs with OK distances and interactions
d = nonzeros(D);         % distances below the large cutoff

% --------- Screen according to maximum difference in nucleotide numbers

if isfield(Query,'MaxDiffMat'),
  if Query.MaxDiffMat(p,q) < Inf,
    k = find(abs(i-j) <= Query.MaxDiffMat(p,q)); %retain pairs close enough together
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to minimum difference in nucleotide numbers

if isfield(Query,'MinDiffMat'),
  if Query.MinDiffMat(p,q) > 1,
    k = find(abs(i-j) >= Query.MinDiffMat(p,q)); %retain pairs far enough apart
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to sign of nucleotide number difference

if isfield(Query,'DifferenceSignMat'),
  if Query.DifferenceSignMat(p,q) < 0,
    k = find(j > i);                          % retain pairs with num2>num1
    i = i(k);
    j = j(k);
    d = d(k);
  elseif Query.DifferenceSignMat(p,q) > 0,
    k = find(j < i);                          % retain pairs with num2<num1
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to the nucleotide mask

if (min(Query.OKCodes{p}) == 0) | (min(Query.OKCodes{q}) == 0),
  if (min(Query.OKCodes{p}) == 1) & (min(Query.OKCodes{q}) == 0),
    k = find(Query.OKCodes{q}(Codes(j)));
  elseif (min(Query.OKCodes{p}) == 0) & (min(Query.OKCodes{q}) == 1),
    k = find(Query.OKCodes{p}(Codes(i)));
  else
    k = find(Query.OKCodes{p}(Codes(i)) .* Query.OKCodes{q}(Codes(j)));
  end

  i = i(k); 
  j = j(k);
  d = d(k);
end

% --------- Screen according to configuration (syn or anti)

if length(Query.Config{p}) > 0,
  switch Query.Config{p}
    case 'syn'
      k = find(cat(1,File.NT(i).Syn) == 1);
    case 'anti'
      k = find(cat(1,File.NT(i).Syn) == 0);
    otherwise
      k = 1:length(i);
  end

  i = i(k);
  j = j(k);
  d = d(k);
end

if length(Query.Config{q}) > 0,
  switch Query.Config{q}
    case 'syn'
      k = find(cat(1,File.NT(j).Syn) == 1);
    case 'anti'
      k = find(cat(1,File.NT(j).Syn) == 0);
    otherwise
      k = 1:length(j);
  end

  i = i(k);
  j = j(k);
  d = d(k);
end

% --------- Screen according to pairwise distance in model

if (Query.Geometric > 0),

  % --------- Compute square of distance difference from model

  d = (d - Query.Distance(p,q)).^2;   % squared difference in distances

  d = d + 0.00000001 * (d == 0);       % avoid rejecting model; make d nonzero

  % --------- Impose upper limit on distance differences; 2-nucleotide cutoff


  if Query.NumNT > 2,
    Wp = Query.LocWeight(p);
    Wq = Query.LocWeight(q);
    MaxD = (Wp + Wq) * (Query.NumNT * Query.DiscCutoff)^2 / (Wp * Wq);
  else
    Wp = 1;
    Wq = 1;
    MaxD = (Query.NumNT * Query.DiscCutoff)^2;
  end

  if isfield(Query,'Flex'),
    MaxD = max(MaxD,Query.Flex(p,q)^2);  % allow larger distance if desired
  end

  k = find(d <= MaxD);            % keep ones with small difference from model

  i = i(k);
  j = j(k);
  d = d(k) * Wp * Wq;

end

% --------- Construct sparse matrix of retained distance difference squares

% [length(i) length(j) length(d) max(k)];

Screen = sparse(i,j,d,File.NumNT,File.NumNT,length(i));



