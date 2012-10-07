% zUpdateExemplars should be used only after zClassLimits has
% been modified.  It redoes the classifications, finds new exemplars, and
% then updates the distance from each pair to the nearest exemplars.

t = cputime;

for f=1:length(File),
  File(f) = zClassifyPairs(File(f));
  File(f).Modified = 0;
  zSaveNTData(File(f));
end

t1 = cputime;

zFindExemplars;

t2 = cputime;

for f=1:length(File),
  File(f) = zUpdateDistanceToExemplars(File(f));
  zSaveNTData(File(f));
end

fprintf('%5.1f minutes to classify pairs\n', (t1-t)/60);
fprintf('%5.1f minutes to find exemplars\n', (t2-t1)/60);
fprintf('%5.1f minutes to update distances to exemplars\n', (cputime-t2)/60);