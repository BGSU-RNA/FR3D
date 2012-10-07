% Direction can be +1 or -1; it tells the order in which to put the 1st cand

function [void] = xAlignCandidates(File,Search,Direction)

Model       = Search.Query;
Candidates  = Search.Candidates;

[L,t] = size(Candidates);

[y,p] = sort(Direction*double(Candidates(1,1:Model.NumNT)));    
                                    % put in increasing or decreasing order
Cand = double(Candidates(:,p));               % re-order nucleotides
F    = Candidates(:,Model.NumNT+1);           % file numbers

if isfield(Model,'MaxDiffMat'),
  MaxDiff = diag(Model.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,Model.NumNT-1);
end

if Direction > 0,
  MaxDiff = MaxDiff;
else
  MaxDiff = fliplr(MaxDiff);
end

maxinsert = zeros(1,Model.NumNT-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(Cand(c,1:Model.NumNT)))-1);
end

% ---------------------------- Print header line

fprintf('               ');
if Model.Geometric > 0,
  fprintf('           ');
end
for j=1:Model.NumNT,
  fprintf('       ');
end
fprintf('    ');
for n = 1:(Model.NumNT-1),
  fprintf('%d',mod(n,10));
  if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
    for i=1:maxinsert(n),
      fprintf(' ');
    end
  else
    fprintf('    ');
  end
end
fprintf('%d\n', mod(Model.NumNT,10));

% ----------------------------- Print alignment

for c = 1:L,                                      % loop through candidates
  f = F(c);                                       % file number
  fprintf('%15s', File(f).Filename);
  if Model.Geometric > 0,
    if isfield(Search,'Discrepancy'),
      fprintf('%11.4f',Search.Discrepancy(c));
    elseif isfield(Search,'AvgDisc'),
      fprintf('%11.4f',Search.AvgDisc(c));
    end
  end
  for j=1:Model.NumNT,                            % print candidate
    fprintf('%3s',File(f).NT(Cand(c,j)).Base);    
    fprintf('%4s',File(f).NT(Cand(c,j)).Number);    
  end
  fprintf('    ');
  for n = 1:(Model.NumNT-1),                      % print alignment
    fprintf('%s', File(F(c)).NT(Cand(c,n)).Base);
    if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
      if Cand(c,n+1) - Cand(c,n) > 1,             % increasing order
        for i = (Cand(c,n)+1):(Cand(c,n+1)-1),
          fprintf('%s', File(F(c)).NT(i).Base);   % show insertions
        end
      elseif Cand(c,n+1) - Cand(c,n) < -1,        % decreasing order
        for i = (Cand(c,n)-1):-1:(Cand(c,n+1)+1),
          fprintf('%s', File(F(c)).NT(i).Base);   % show insertions
        end
      end
      for i=1:(1 + maxinsert(n) - abs(Cand(c,n+1)-Cand(c,n))),
        fprintf('-');
      end
    else
      fprintf('....');
    end
  end
  fprintf('%s', File(F(c)).NT(Cand(c,Model.NumNT)).Base);
  fprintf('\n');
  drawnow
end

