
function [Screen] = xPairwiseScreen(File,Codes,Model,p,q,PC);

% --------- Screen according to interaction between nucleotides

if Model.Geometric == 0,
  D = File.Distance .* (File.Distance < Model.Diameter);
                                        % cap distance for non-geometric search
else
  D = File.Distance;
end

if isfield(Model,'EdgeNums'),                 % if screening by edges
  if length(Model.EdgeNums{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Model.EdgeNums{p,q}),
      E = E + (fix(File.Edge) == Model.EdgeNums{p,q}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Model,'ExcludeEdges'),                 % if screening by edges
  if length(Model.ExcludeEdges{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Model.ExcludeEdges{p,q}),
      E = E + (fix(File.Edge) == Model.ExcludeEdges{p,q}(i));
    end
    D = D .* (E == 0);
  end
end
    
if isfield(Model,'OKPairs'),                 % if screening by paircode
  if length(Model.OKPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Model.OKPairs{p,q}),
      E = E + sparse(PC == Model.OKPairs{p,q}(i));
    end
    D = D .* (E > 0);
  end
end
    
if isfield(Model,'ExPairs'),                 % if screening by paircode
  if length(Model.ExPairs{p,q} > 0),
    E = sparse(zeros(size(D)));
    for i=1:length(Model.ExPairs{p,q}),
      E = E + sparse(PC == Model.ExPairs{p,q}(i));
    end
    D = D .* (E == 0);
  end
end
    
[i,j] = find(D);         % nucleotide pairs with OK distances and interactions
d = nonzeros(D);         % distances below the large cutoff

% --------- Screen according to maximum difference in nucleotide numbers

if isfield(Model,'MaxDiffMat'),
  if Model.MaxDiffMat(p,q) < Inf,
    k = find(abs(i-j) <= Model.MaxDiffMat(p,q)); %retain pairs close enough together
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to minimum difference in nucleotide numbers

if isfield(Model,'MinDiffMat'),
  if Model.MinDiffMat(p,q) > 1,
    k = find(abs(i-j) >= Model.MinDiffMat(p,q)); %retain pairs far enough apart
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to sign of nucleotide number difference

if isfield(Model,'DifferenceSignMat'),
  if Model.DifferenceSignMat(p,q) < 0,
    k = find(j > i);                          % retain pairs with num2>num1
    i = i(k);
    j = j(k);
    d = d(k);
  elseif Model.DifferenceSignMat(p,q) > 0,
    k = find(j < i);                          % retain pairs with num2<num1
    i = i(k);
    j = j(k);
    d = d(k);
  end
end

% --------- Screen according to the nucleotide mask

if (min(Model.OKCodes{p}) == 0) | (min(Model.OKCodes{q}) == 0),
  if (min(Model.OKCodes{p}) == 1) & (min(Model.OKCodes{q}) == 0),
    k = find(Model.OKCodes{q}(Codes(j)));
  elseif (min(Model.OKCodes{p}) == 0) & (min(Model.OKCodes{q}) == 1),
    k = find(Model.OKCodes{p}(Codes(i)));
  else
    k = find(Model.OKCodes{p}(Codes(i)) .* Model.OKCodes{q}(Codes(j)));
  end

  i = i(k); 
  j = j(k);
  d = d(k);
end

% --------- Screen according to configuration (syn or anti)

if length(Model.Config{p}) > 0,
  switch Model.Config{p}
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

if length(Model.Config{q}) > 0,
  switch Model.Config{q}
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

if (Model.Geometric > 0),

  % --------- Compute square of distance difference from model

  d = (d - Model.Distance(p,q)).^2;   % squared difference in distances

  d = d + 0.00000001 * (d == 0);       % avoid rejecting model; make d nonzero

  % --------- Impose upper limit on distance differences; 2-nucleotide cutoff


  if Model.NumNT > 2,
    Wp = Model.LocWeight(p);
    Wq = Model.LocWeight(q);
    MaxD = (Wp + Wq) * (Model.NumNT * Model.DiscCutoff)^2 / (Wp * Wq);
  else
    Wp = 1;
    Wq = 1;
    MaxD = (Model.NumNT * Model.DiscCutoff)^2;
  end

  if isfield(Model,'Flex'),
    MaxD = max(MaxD,Model.Flex(p,q)^2);  % allow larger distance if desired
  end

  k = find(d <= MaxD);            % keep ones with small difference from model

  i = i(k);
  j = j(k);
  d = d(k) * Wp * Wq;

end

% --------- Construct sparse matrix of retained distance difference squares

% [length(i) length(j) length(d) max(k)];

Screen = sparse(i,j,d,File.NumNT,File.NumNT,length(i));



