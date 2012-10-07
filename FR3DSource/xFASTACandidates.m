% xFASTACandidates(File,Search,Direction) shows a multiple sequence
% alignment of the candidates in Search.  Bases which correspond in the
% geometric search are aligned with one another.  If a maximum distance
% has been specified, or the maximum gap is small, bases between the
% bases in the candidate are also displayed.

% Direction can be +1 or -1; it tells the order in which to put the
% nucleotides in the first candidate

function [Text] = xFASTACandidates(File,Search,Direction,Param)

Query      = Search.Query;
Candidates = Search.Candidates;
N          = Query.NumNT;

[L,t] = size(Candidates);

[y,p] = sort(Direction*double(Candidates(1,1:N)));    
                                    % put nucleotides in inc/decreasing order
Cand = double(Candidates(:,p));     % re-order nucleotides
F    = Candidates(:,N+1);           % file numbers

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(p,p),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

if Direction > 0,
  MaxDiff = MaxDiff;
else
  MaxDiff = fliplr(MaxDiff);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
end

% ---------------------------- Print header line

t = 1;                                          % line of text we're on
k = N;                                          % number of columns being shown
Text{t} = sprintf('               ');
if Query.Geometric > 0,
  Text{t} = [Text{t} sprintf('           ')];
end
for j=1:N,
  Text{t} = [Text{t} sprintf('        ')];
end
Text{t} = [Text{t} sprintf('    ')];
for n = 1:(N-1),
  Text{t} = [Text{t} sprintf('%d',mod(n,10))];
  if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
    for i=1:maxinsert(n),
      Text{t} = [Text{t} sprintf(' ')];
      k = k + 1;                                % more columns being shown
    end
  else
    Text{t} = [Text{t} sprintf('    ')];
  end
end
Text{t} = [Text{t} sprintf('%d', mod(N,10))];

% ----------------------------- Print alignment

CorrCodeList = zeros(c,N);                        % for bases corresp to query
for j = 1:N,
  InsCode{j}.Count = zeros(1,4);
  InsLength{j}.Count = [];
end

for c = 1:L,                                      % loop through candidates
  f = F(c);                                       % file number
  t = t + 1;
  Text{t} = [sprintf('> %s %5s %s%5s %s%5s', Search.Query.Name, File(f).Filename, File(f).NT(Cand(c,1)).Base, File(f).NT(Cand(c,1)).Number, File(f).NT(Cand(c,end)).Base, File(f).NT(Cand(c,end)).Number)];
  t = t + 1;
  Text{t} = '';

  k = 1;                                          % column counter

  for n = 1:(N-1),                                % print alignment for cand
    Text{t} = [Text{t} sprintf('%s', File(F(c)).NT(Cand(c,n)).Base)];
    j = File(F(c)).NT(Cand(c,n)).Code;
    CorrCodeList(c,n) = j;                        % store code of corresp base
    k = k + 1;
    if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
      if Cand(c,n+1) - Cand(c,n) > 1,             % increasing order
        for i = (Cand(c,n)+1):(Cand(c,n+1)-1),
          Text{t} = [Text{t} sprintf('%c', File(F(c)).NT(i).Base)];   % show insertions
          j = File(F(c)).NT(i).Code;
          InsCode{n}.Count(j) = InsCode{n}.Count(j) + 1;
        end
      elseif Cand(c,n+1) - Cand(c,n) < -1,        % decreasing order
        for i = (Cand(c,n)-1):-1:(Cand(c,n+1)+1),
          Text{t} = [Text{t} sprintf('%c', File(F(c)).NT(i).Base)];   % show insertions
          j = File(F(c)).NT(i).Code;
          InsCode{n}.Count(j) = InsCode{n}.Count(j) + 1;
        end
      end
      for i=1:(1 + maxinsert(n) - abs(Cand(c,n+1)-Cand(c,n))),
        Text{t} = [Text{t} sprintf('-')];
        k = k + 1;
      end
    else
      Text{t} = [Text{t} sprintf('*')];
    end

    h = abs(Cand(c,n+1)-Cand(c,n));                % number of insertions + 1
    if h > length(InsLength{n}.Count),
      InsLength{n}.Count(h) = 1;
    else
      InsLength{n}.Count(h) = InsLength{n}.Count(h) + 1;
    end

  end
  Text{t} = [Text{t} sprintf('%s', File(F(c)).NT(Cand(c,N)).Base)];
  j = File(F(c)).NT(Cand(c,N)).Code;
  CorrCodeList(c,N) = j;                           % store last corresp base

  drawnow
end

if nargin < 4,
  for t = 1:length(Text),
    fprintf('%s\n', Text{t});
  end
end
