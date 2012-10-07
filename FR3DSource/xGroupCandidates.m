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
  Lab{i} = [b n ' ' File(f).Filename];
end

D = Search.Disc(Done,Done);                    % mutual distances to consider

[s,t] = size(D);
for i=1:s,
  D(i,i) = 0;                                  % just in case
end

Y = squareform(full(D));                       % convert to a vector
Z = linkage(Y,'average');                      % compute cluster tree

L = length(Search.Candidates(:,1));
ClusterNum = L+1;
CurrentCluster = 1:L;

for i = 1:L,
  Label{i} = '111                   ';
end

for pp = 1:length(Z(:,1));
  i = find(CurrentCluster == Z(pp,1));
  for ii = 1:length(i),
    Label{i(ii)} = ['1' Label{i(ii)}];
    CurrentCluster(i(ii)) = ClusterNum;
  end

  j = find(CurrentCluster == Z(pp,2));
  for jj = 1:length(j),
    Label{j(jj)} = ['2' Label{j(jj)}];
    CurrentCluster(j(jj)) = ClusterNum;
  end

  ClusterNum = ClusterNum + 1;
end

for i = 1:L,
  La = Label{i};
  SLabel{i} = '';                       % short label
  for j = 1:fix((length(La)-1)/4),
    S = La((1+4*(j-1)):(4*j));
    if any(S == ' '),
      T = La((1+4*(j-1)):end);
    else
      switch S,
      case '1111', T = 'A';
      case '1112', T = 'B';
      case '1121', T = 'C';
      case '1122', T = 'D';
      case '1211', T = 'F';
      case '1212', T = 'G';
      case '1221', T = 'I';
      case '1222', T = 'J';
      case '2111', T = 'Z';
      case '2112', T = 'Y';
      case '2121', T = 'W';
      case '2122', T = 'V';
      case '2211', T = 'T';
      case '2212', T = 'S';
      case '2221', T = 'Q';
      case '2222', T = 'P';
      otherwise, T = 'Z';
      end
    end
    SLabel{i}= [SLabel{i} T];
  end
  FLabel{i} = [Lab{i} ' ' SLabel{i}];
end

figure(99)
[H,T,p] = dendrogram(Z,0,'colorthreshold',0.2,'orientation','left','labels',FLabel);
%set(gca,'FontSize',1)             % helpful when printing to PDF
%set(gcf,'Renderer','painters');

Search.GroupLabel = SLabel;

n = length(p):-1:1;
r = Done(p(n));
q = Done(p);                                      % reverse order

% ----------------------------------- List and display results

fprintf('Candidates are listed in order of groups:\n');
fprintf('\n');

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(q,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(q);
S.Disc         = Search.Disc(q,q);
S.DiscComputed = Search.Disc(1,q);
S.File         = Search.File;
S.GroupLabel   = Search.GroupLabel(q);
if isfield(Search,'AvgDisc'),
  S.AvgDisc = Search.AvgDisc(q);
end

S.CandidateFilenames = Search.CandidateFilenames;

xListCandidates(S,Inf);                             % show on screen

xDisplayCandidates(File,S,Level+1,UsingFull);       % display, level 1
