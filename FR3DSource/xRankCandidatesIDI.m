% xRankCandidatesIDI(File,Model,Candidates) calculates the discrepancy for
% each candidate; it is only for two nucleotide candidates right now!

function [Discrepancy, Candidates, i] = xRankCandidatesIDI(File,Model,Cand,Verbose)

if nargin < 4,
  Verbose = 0;
end

if Verbose > 0,
  fprintf('Calculating discrepancy\n');
end

N = length(Model.NT);
s = length(Cand(:,1));

Discrepancy = zeros(s,1);
Candidates  = uint16(zeros(s,length(Cand(1,:))));

count = 0;

tic

if Verbose > 0,
  fprintf('Seconds remaining:');
end


for i=1:s,
%  A = xDiscrepancyFast(Model,File(Cand(i,N+1)).NT(Cand(i,[1:N])));
  f = Cand(i,N+1);
  A = zIsoDiscrepancy(Model.NT(1),Model.NT(2),File(f).NT(Cand(i,1)),File(f).NT(Cand(i,2)));

  if A >= 0,
    count = count + 1;
    Discrepancy(count,1) = A;
    Candidates(count,:)  = Cand(i,:);
  end

  if (mod(i,round(s/10)) == 0) && (Verbose > 0)
    fprintf(' %d', fix((s-i)*toc/i)); 
    drawnow
  end
end

if Verbose > 0,
  fprintf('\n');
end

Candidates  = Candidates(1:count,:);

[y,i]       = sort(Discrepancy);                    % sort by discrepancy
Candidates  = Candidates(i,:);
Discrepancy = Discrepancy(i);

if Verbose > 1,
  fprintf('Calculating discrepancy took        %8.3f seconds\n',toc);
end

drawnow

