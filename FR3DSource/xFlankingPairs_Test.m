
Filenames = {'1J5E','2AVY','1S72','3U5H','1M5P','1M5V','2E5L','2GJW','357D','364D','3J0P','3J0Q','3OL6','3OL7','3OL8','3OL9','3OLA','3OLB'};

for f = 1:length(Filenames),

File = zAddNTData(Filenames{f});

FF = xFlankingPairs(File);
TT = triu(FF.Flank);

F = xFlankingPairs_new(File); 
T = triu(F.Flank);


if sum(sum(T > TT)) > 0,
  fprintf('Some new flanking pairs were introduced!!!\n');
end


if sum(sum(TT > T)) > 0,
  fprintf('Some flanking pairs have been dropped from %s.  They are:\n', File.Filename);
  [i,j] = find(TT > T);
  for k = 1:length(i),
    N1 = File.NT(i(k));
    N2 = File.NT(j(k));
    fprintf('%s%s_%s with %s%s_%s   \n', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain);
  end
end

if sum(sum(abs(TT-T))) == 0,
  fprintf('No change in flanking pairs in %s.\n', File.Filename);
end

fprintf('\n');

end
