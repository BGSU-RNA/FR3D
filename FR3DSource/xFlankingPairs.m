% xFlankingPairs identifies nucleotides on the same strand which make 
function [File] = xFlankingPairs(File)

E = triu(fix(abs(File.Edge))==1) .* (File.Range == 0); % nested cWW's, i < j
X = (File.Range == 0) .* (abs(File.Edge) > 0) .* (abs(File.Edge) < 15);

H = sparse(zeros(File.NumNT,File.NumNT));              % to indicate flanks

[i,j] = find(E);                      % indices of NT's making nested cWW's

a = [i; j];                           % all indices of nested cWW pairs

a = sort(a);                          % put them into one list

for k = 1:(length(a)-1),              % run through all NT making nested cWW
  if a(k+1) - a(k) > 1,               % where there is a gap in the number,
    H(a(k),a(k+1)) = 1;               % these two flank something
  end
end

H = H + H';                           % symmetrize H

File.Flank = H;                       % store

