% xReduceRedundantCandidates(Candidates,Discrepancy) removes candidates which are redundant with each other because they occur in redundant chains

% load 2009-07-18_14_33_04-Redundant_chains_geometric
% Candidates = Search.Candidates;
% Discrepancy = Search.Discrepancy;

function [Candidates,Discrepancy] = xReduceRedundantCandidates(Candidates,Discrepancy)

N  = length(Candidates(1,:)) - 1;  % number of nucleotides
OK = ones(size(Candidates(:,1)));  % rows of Candidates which are kept

Filenums = unique(Candidates(:,N+1));  % file numbers which occur

M = sparse(eye(max(max(Candidates(:,1:N)))));

for i = 1:length(Filenums),
  f = Filenums(i);                     % current file number
  c = find(Candidates(:,N+1) == f);    % candidates from the same file
  S = sparse(zeros(length(c),length(c)));
  for n = 1:N,
    S = S + M(Candidates(c,n),Candidates(c,n));
  end

  for cc = 1:length(c),
    S(cc,cc) = 0;                      % cand is not redundant with itself
  end

  for cc = 1:length(c),
    d = find(S(cc,:) > N/2);
    OK(c(d)) = zeros(1,length(d));
    S(d,:) = zeros(length(d),length(c));
  end
end

Candidates  = Candidates(find(OK),:);
Discrepancy = Discrepancy(find(OK));
