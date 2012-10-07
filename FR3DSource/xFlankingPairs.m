% xFlankingPairs identifies nucleotides which flank a single-stranded region by making nested canonical cWW basepairs

function [File] = xFlankingPairs(File)

if File.NumNT > 0,
  H = sparse([],[],[],File.NumNT,File.NumNT);        % matrix to indicate flanks

  E = triu(fix(abs(File.Edge))==1);     % cWW pairs
  [i,j] = find(E);                      % indices of NT's making cWWs

  % Paircode list
  % 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
  % 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

  %  paircode = 4*(N2.Code-1) + N1.Code;           % AA is 1, CA is 2, etc.

  Nested = zeros(1,length(i));          % place to stored nested canonical cWWs
  for u = 1:length(i),
    pc = 4*(File.NT(j(u)).Code-1) + File.NT(i(u)).Code;
    if File.Crossing(i(u),j(u)) == 0 && any(pc == [4 7 10 12 13 15]), % canonical
      Nested(u) = 1;
    end
  end

  k = find(Nested);
  i = i(k);
  j = j(k);

  a = [i; j];                           % all indices of nested canonical cWWs

  a = sort(a);                          % put them into one list

  for k = 1:(length(a)-1),              % run through all NT making nested cWW
    if a(k+1) - a(k) > 1,               % where there is a gap in the number,
      H(a(k),a(k+1)) = 1;               % these two flank a single-stranded region
    end
  end

  H = H + H';                           % symmetrize H
 
  File.Flank = H;                       % store

else

  File.Flank = [];

end


if 0 > 1,
  D = (File.Flank - F.Flank);
  [i,j] = find(D == 1);               % was flank, not now
  for k = 1:length(i),
    b = find(fix(abs(File.Edge(i(k),:))) == 1);
    fprintf('cWW pair was %s%4s-%s%4s\n', File.NT(i(k)).Base, File.NT(i(k)).Number, File.NT(b(1)).Base, File.NT(b(1)).Number);
  end

  [i,j] = find(D == -1);               % was flank, not now
  for k = 1:length(i),
    b = find(fix(abs(File.Edge(i(k),:))) == 1);
    fprintf('New cWW pair is %s%4s-%s%4s\n', File.NT(i(k)).Base, File.NT(i(k)).Number, File.NT(b(1)).Base, File.NT(b(1)).Number);
  end
end

