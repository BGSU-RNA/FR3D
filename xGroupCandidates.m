% xGroupCandidates clusters sequences and displays groups together

function [Search] = xGroupCandidates(File,Search,Level,UsingFull)

N      = Search.Query.NumNT;
Search = xMutualDiscrepancy(File,Search);      % compute discrepancy matrix
Done   = find(Search.DiscComputed);            % ones already computed

for i=1:length(Done),
  f = Search.Candidates(Done(i),N+1);          % file number
  b = File(f).NT(Search.Candidates(Done(i),1)).Base;
  n = File(f).NT(Search.Candidates(Done(i),1)).Number;
  n = sprintf('%4s',n);
  if Search.Query.Geometric > 0,
      if isfield(Search,'AvgDisc'),
        d = sprintf('%6.4f',Search.AvgDisc(Done(i)));
      else
        d = sprintf('%6.4f',Search.Discrepancy((Done(i))));
      end
    else
      d = sprintf('%5d',Search.Discrepancy(Done(i))); % orig candidate number
    end
  Lab{i} = [b n ' ' num2str(d)];
end

D = Search.Disc(Done,Done);                    % mutual distances to consider

[s,t] = size(D);
for i=1:s,
  D(i,i) = 0;                                  % just in case
end

Y = squareform(full(D));                       % convert to a vector
Z = linkage(Y,'average');                      % compute cluster tree
figure(25)
[H,T,p] = dendrogram(Z,0,'colorthreshold',0.2,'orientation','left','labels',Lab);

n = length(p):-1:1;
r = Done(p(n));
q = Done(p);                                      % reverse order

% ----------------------------------- List and display results

fprintf('Candidates sorted by centrality within these candidates:\n');
fprintf('\n');

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(r,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(r);
S.Disc         = Search.Disc(r,r);
S.DiscComputed = Search.Disc(1,r);
S.File         = Search.File;
S.CandidateFilenames = Search.CandidateFilenames;
if isfield(Search,'AvgDisc'),
  S.AvgDisc = Search.AvgDisc(r);
end

xListCandidates(S,Inf);                 % show on screen

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(q,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(q);
S.Disc         = Search.Disc(q,q);
S.DiscComputed = Search.Disc(1,q);
S.File         = Search.File;
if isfield(Search,'AvgDisc'),
  S.AvgDisc = Search.AvgDisc(q);
end

xDisplayCandidates(File,S,Level+1,UsingFull);          % display, level 1
