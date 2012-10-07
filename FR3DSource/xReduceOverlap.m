% xReduceOverlap(Candidates,Discrepancy,Verbose) 

function [Candidates, Discrepancy] = xReduceOverlap(Candidates, Discrepancy, Verbose)

if nargin < 3,
  Verbose = 0;
end

N = length(Candidates(1,:)) - 1;                 % number of nucleotides

if Verbose > 0,
  fprintf('%6d candidates after %2d exclusions\n', length(Candidates(:,1)),0);
  drawnow
end

a = 1;

for a=1:N,
  i = [a:N 1:(a-1)];                             % cyclic permutations
  [Candidates, Discrepancy] = xExcludeOverlap2(Candidates, Discrepancy, i);
  if Verbose > 0,
    fprintf('%6d candidates after %2d exclusions\n',length(Candidates(:,1)),a);
    drawnow
  end
end

for a=1:N,
  i = [a:-1:1 N:-1:(a+1)];                       % more permutations
  [Candidates, Discrepancy] = xExcludeOverlap2(Candidates, Discrepancy, i);
  %fprintf('%6d candidates after %2d exclusions\n',length(Candidates(:,1)),a+N);
  drawnow
end

return

while a < 3000,
  r = rand(1,N);
  [y,i] = sort(r);                    % i is a random permutation
    
  [Candidates, Discrepancy] = xExcludeOverlap2(Candidates, Discrepancy, i);

  [a size(Candidates)]
  drawnow

  a = a + 1;
end
