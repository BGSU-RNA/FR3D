% xNeighborhood(File,Indices,v,MaxDiff,MaxInsert) returns indices "near"
% the given indices, in a way determined by the code v

function [NewIndices] = xNeighborhood(File,Indices,v,Display)

MaxDiff = Display.MaxDiff;
MaxInsert = Display.MaxInsert;

p = Display.p;

N = length(Indices);

switch v,
  case 0,
    NewIndices = Indices;
  case 1,
    NewIndices = Indices;
    for n = 1:(N-1),
      if (MaxDiff(n) < Inf) | (MaxInsert(n) < 16),   % if only few insertions
        NewIndices = [NewIndices (Indices(n)+1):(Indices(n+1)-1)];
      end
    end
  case 2,
    NewIndices = Indices;
    for n = 1:(N-1),
      if (MaxDiff(n) < Inf) | (MaxInsert(n) < 20),   % if only few insertions
        NewIndices = [NewIndices (Indices(n)+1):(Indices(n+1)-1)];
      end
    end
  otherwise,
    d = [1 1 6 8 10];
    a = zeros(1,File.NumNT);
    for j=1:length(Indices),
      a = a + (File.Distance(Indices(j),:) < d(v)) .* ...
              (File.Distance(Indices(j),:) > 0);
    end
    a(Indices) = zeros(1,length(Indices));  % take out ones in Indices
    B = find(a);
    NewIndices = [Indices B];
end

NewIndices = sort(NewIndices);

% it should make sure there are no duplicates in NewIndices
% it should check File to make sure that everything in NewIndices really
% is present in File, otherwise something will crash
