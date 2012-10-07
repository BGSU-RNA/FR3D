
% histogram distance to nearest exemplar for classified pairs

d = [];
c = 1;

for f = 1:length(File),
 for p = 1:length(File(f).Pair),
  Pair = File(f).Pair(p);
  if (abs(Pair.Edge) < 20) && (abs(Pair.Edge) > 0),
   if fix(abs(Pair.Edge)) ~= fix(abs(Pair.Classes(1))),
%   if (abs(Pair.Edge)) ~= (abs(Pair.Classes(1))),
    d(c,:) = [Pair.Distances(1) Pair.Paircode Pair.Edge Pair.Classes(1)];
    c = c + 1;

    fprintf('File %s Nucleotides %s%4s and %s%4s Classified as %4s %8.2f but  closest to %4s %8.2f\n', File(f).Filename, File(f).NT(Pair.Base1Index).Base, File(f).NT(Pair.Base1Index).Number, File(f).NT(Pair.Base2Index).Base, File(f).NT(Pair.Base2Index).Number, zEdgeText(Pair.Edge), Pair.Edge, zEdgeText(Pair.Classes(1)), Pair.Classes(1));
%Pair
   end
  end
 end
end

clf
hist(d(:,1),30);
title('Distance to nearest exemplar');

dd = sortrows(d, 3);
