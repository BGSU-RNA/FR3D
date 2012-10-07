% zUpdateDistanceToExemplars computes the distance between each pair in
% File and the exemplars for each category of interaction.  This indicates
% how close the pair is to the top three categories.

function [File] = zUpdateDistanceToExemplars(File)

if exist('PairExemplars.mat','file') > 0,

  load('PairExemplars','Exemplar');

  for f = 1:length(File),
    for p = 1:length(File(f).Pair),
      Pair = File(f).Pair(p);
      NT1 = File(f).NT(Pair.Base1Index);
      NT2 = File(f).NT(Pair.Base2Index);
      [c,d,h] = zDistanceToExemplars(Exemplar,NT1,NT2);
      File(f).Pair(p).Classes   = c(1:3);
      File(f).Pair(p).Distances = d(1:3);
    end
%    fprintf('Updated %10s\n', File(f).Filename);
%    zSaveNTData(File(f));
  end
else
  for p = 1:length(File.Pair),                  % fill in fictitious distances
    File.Pair(p).Classes   = 99 * ones(1,3);
    File.Pair(p).Distances = 99999999 * ones(1,3);
    File.Pair(p).ExemIndex = [1 2 3];
  end
end
