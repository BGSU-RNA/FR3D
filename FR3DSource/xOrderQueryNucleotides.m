%  zOrderQueryNucleotides finds a reasonably good order to put the query nucleotides to speed up xFindPolyhedra

function [Perm] = xFindPermutation(NNZ);

[m,n] = size(NNZ);

if m <= 3,
  Perm = 1:m;
else

  minp = Inf;
  mini = 1;
  minj = 2;
  mink = 3;

  for i = 1:n,
    for j = (i+1):n,
      for k = (j+1):n,
        p = NNZ(i,j) * NNZ(i,k) * NNZ(j,k);
        if p < minp,
          minp = p;
          mini = i;
          minj = j;
          mink = k;
        end
      end
    end
  end

  i = mini;
  j = minj;
  k = mink;

  a = [NNZ(i,j) NNZ(i,k) NNZ(j,k)];

  [y,z] = min(a);

  Chosen = [i j k];

  switch z(1)
  case 1
    Chosen = [i j k];
  case 2
    Chosen = [i k j];
  case 3
    Chosen = [j k i];
  end

  Remain = setdiff(1:n,Chosen);

  while ~isempty(Remain),
    minp = Inf;
    minr = 1;
    for r = 1:length(Remain)
      p = prod(NNZ(Chosen,Remain(r)));           % product 
      if p < minp,
        minp = p;
        minr = r;
      end
    end

    Chosen = [Chosen Remain(minr)];
    Remain = setdiff(Remain,Remain(minr));
  end

  Perm = Chosen;
end
