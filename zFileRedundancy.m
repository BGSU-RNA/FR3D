% zFileRedundancy(File) explores possible redundancy between PDB files

function [void] = zFileRedundancy(File)

a = zeros(length(File),7);

for f=1:length(File),
  a(f,1) = length(File(f).NT);
  a(f,2) = length(File(f).Pair);
end

[A,i] = sortrows(a,[-1 -2]);


B = [];

for j = 1:length(i),
 fprintf('%5s has %4d nucleotides, %5d pairs, ', File(i(j)).Filename, A(j,1), A(j,2));
 if isempty(File(i(j)).Info.Resolution),
   fprintf('resolution  ????');
 else
   fprintf('resolution %5.2f', File(i(j)).Info.Resolution);
 end

 if A(j,1) > 0,

  B = [];
  for k = 1:A(j,1),
    B = [B File(i(j)).NT(k).Base];
  end
  
  if j > 1,
    if length(B) == length(PreviousB),
      d = xDiscrepancy(File(i(j)),1:A(j,1),File(i(j-1)),1:A(j-1,1));
      fprintf(' %7.4f geometric discrepancy',d);

      p = sum(B == PreviousB)/length(B);
      if p > 0.9,
        fprintf(',%4.0f%% base agreement with previous file.', 100*p);
      else
        fprintf(' with previous file                      ');
      end
    else
      fprintf('                                                                       ');
    end
  else
    fprintf('                                                                       ');
  end

  Info = File(i(j)).Info;

  fprintf('%20s %40s %40s %70s', Info.Type, Info.RNA, Info.Species, Info.LigandsAndComments);

  fprintf('\n');

  PreviousB = B;
 else
  fprintf('\n');
 end
end


