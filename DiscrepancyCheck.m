
% histogram distance to nearest exemplar for classified pairs

d = [];
c = 1;

for f = 1:length(File),
 for p = 1:length(File(f).Pair),
  Pair = File(f).Pair(p);
%  if (abs(Pair.Edge) < 20) && (abs(Pair.Edge) > 0),
  if (abs(Pair.Edge) > 0),
    d(c,:) = [Pair.Distances(1) Pair.Paircode Pair.Edge Pair.Classes(1)];
    c = c + 1;
  end
 end
end

clf
hist(d(:,1),30);
title('Distance to nearest exemplar');

dd = dd(find(d(:,3) == 30))
dd = sortrows(d, 3);