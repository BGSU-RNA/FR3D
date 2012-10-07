% zCountPairs counts the number of nucleotides involved in different types of interactions

%cWW only:   nucleotides making a cWW pair but no other basepair

%cWW and 1 ncWW: nucleotides making both a cWW pair and exactly one non-cWW pair

%cWW and 2+ ncWW: nucleotides making both a cWW pair and two or more non-cWW pairs

%1 ncWW: nucleotides making no cWW pair, but exactly one ncWW pair

%2 ncWW: nucleotides making no cWW pair, but exactly two ncWW pairs

%3 ncWW: nucleotides making no cWW pair, but exactly three ncWW pairs

%4+ ncWW: nucleotides making no cWW pair, but four or more ncWW pairs

%stack only: nucleotides making no basepairs, but at least one stacking

%nothing: nucleotides making no basepair and no stacking pair

%Check total number: compares the total of the above with the known number of nucleotides in the structure, to make sure nothing is missing


% File = zAddNTData('1s72',2);             % load file data in this way

function [void] = zCountPairs(File,Multiple)

if nargin < 2,
  Multiple = 0;
end

fprintf('FileID cWW only - cWW and 1 ncWW - cWW and 2+ ncWW - 1 ncWW - 2 ncWW - 3 ncWW - 4+ ncWW - stack only - nothing Check total number\n');

for f = 1:length(File),
  A = 0;
  B = 0;
  C = 0;
  D = zeros(1,10);
  G = 0;
  H = 0;
  E = fix(abs(File(f).Edge));
  cWW   = sum(E == 1, 2);
  ncWW  = sum((1 < E) .* (E < 15),2);
  stack = sum((19 < E) .* (E < 24),2);
  for n = 1:length(File(f).NT),
    if (cWW(n) > 0) && (ncWW(n) == 0),
      A = A + 1;
    elseif (cWW(n) > 0) && (ncWW(n) == 1),
      B = B + 1;
    elseif (cWW(n) > 0) && (ncWW(n) > 1),
      H = H + 1;

      if (ncWW(n) > 1) && (Multiple > 0),
        fprintf('%s nucleotide %s%s makes these pairs: ', File(f).Filename, File(f).NT(n).Base, File(f).NT(n).Number);
        P = find((1 <= E(n,:)) .* (E(n,:) < 15));
        for i = 1:length(P),
          fprintf('%s with %s%s, ', zEdgeText(File(f).Edge(n,P(i))), File(f).NT(P(i)).Base,File(f).NT(P(i)).Number);
        end
        fprintf('\n');
      end

    elseif (cWW(n) == 0) && (ncWW(n) == 0) && (stack(n) > 0),
      C = C + 1;
    elseif (cWW(n) == 0) && (ncWW(n) == 0) && (stack(n) == 0),
      G = G + 1;
    elseif (cWW(n) == 0) && (ncWW(n) > 0),
      D(ncWW(n)) = D(ncWW(n)) + 1;

      if (ncWW(n) > 2) && (Multiple > 0),
        fprintf('%s nucleotide %s%s makes these pairs: ', File(f).Filename, File(f).NT(n).Base, File(f).NT(n).Number);
        P = find((1 < E(n,:)) .* (E(n,:) < 15));
        for i = 1:length(P),
          fprintf('%s with %s%s, ', zEdgeText(File(f).Edge(n,P(i))), File(f).NT(P(i)).Base,File(f).NT(P(i)).Number);
        end
        fprintf('\n');
      end

    end
  end
  fprintf('%s   %8d         %8d          %8d %8d %8d %8d  %8d     %8d  %8d %8d %8d\n', File(f).Filename, A, B, H, D(1), D(2), D(3), sum(D(4:end)), C, G, A+B+H+sum(D)+C+G, length(File(f).NT));

end

