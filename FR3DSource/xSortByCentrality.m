% xSortByCentrality re-orders the candidates according to their centrality

function [Search] = xSortByCentrality(File,Search,Level,UsingFull)

Search = xMutualDiscrepancy(File,Search);

s = length(find(Search.DiscComputed));            % number computed

% ----------------------------------- Sort by centrality

%[z,j] = sort(sum(Search.Disc));               % sort by average discrepancy
[z,j] = sort(max(Search.Disc));              % sort by maximum discrepancy

z = z / (s-1);                                % average discrepancy among these

% ----------------------------------- List and display results

fprintf('Candidates sorted by centrality within these candidates:\n');
fprintf('\n');

S.Query        = Search.Query;                      % set up new "Search" data
S.Candidates   = Search.Candidates(j,:);            % re-order candidates
S.Discrepancy  = Search.Discrepancy(j);
S.Disc         = Search.Disc(j,j);
S.DiscComputed = Search.Disc(1,j);
S.AvgDisc      = z;
S.DisttoCenter = Search.Disc(j(1),j);
S.File         = Search.File;
S.CandidateFilenames = Search.CandidateFilenames;

xListCandidates(S,Inf);                             % show on screen
xDisplayCandidates(File,S,Level+1,UsingFull);       % display, level 1

