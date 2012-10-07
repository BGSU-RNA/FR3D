% xRankCandidates(File,Model,Candidates) calculates the discrepancy for
% each candidate

function [Discrepancy, Candidates] = xRankCandidates(File,Model,Cand,Verbose)

if nargin < 4,
  Verbose = 0;
end

if Verbose > 0,
  fprintf('Calculating discrepancy\n');
end

N = Model.NumNT;
s = length(Cand(:,1));

Discrepancy = zeros(65000,1);
Candidates  = uint16(zeros(65000,length(Cand(1,:))));

count = 0;

tic

if Verbose > 0,
  fprintf('Seconds remaining:');
end


for i=1:s,
  A = xDiscrepancyFast(Model,File(Cand(i,N+1)).NT(Cand(i,[1:N])));

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
Discrepancy = sqrt(Discrepancy(1:count,1))/Model.NumNT;

[y,i]       = sort(Discrepancy);                    % sort by discrepancy
Candidates  = Candidates(i,:);
Discrepancy = Discrepancy(i);

if Verbose > 1,
  fprintf('Calculating discrepancy took        %8.3f seconds\n',toc);
end

drawnow

