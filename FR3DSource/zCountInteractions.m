% zCountInteractions(File) tabulates the number of each type of interaction in File
% It formats the counts by geometric family and, within each family, by basepair.
% It writes an Excel sheet with a standard format.
% CountFilename could be 'Basepair_counts_NonRedundant_2008_02_21_list.xls'
% If omitted, the filenames in File are used to construct the filename.

function [void] = zCountInteractions(File,CountFilename)

if nargin == 1,
  CountFilename = 'Basepair_counts';
  for f = 1:length(File),
    CountFilename = [CountFilename '_' File(f).Filename];
  end
  CountFilename = [CountFilename '.xls'];
end

MaxCat = 25;

CI = [];

for f = 1:length(File),
  E = File(f).Edge;
  E = E .* (E > 0) .* (E < MaxCat);           % only certain interactions
  [i,j,k] = find(E);
  Code1 = cat(1,File(f).NT(i).Code);
  Code2 = cat(1,File(f).NT(j).Code);
  
  CI = [CI; [k Code1 Code2]];
end

CI = [floor(CI(:,1)*10) CI(:,2:3)];           % multiply classes by 10

Count = zeros(4,4,MaxCat*10);

for i = 1:length(CI(:,1)),
  Count(CI(i,2),CI(i,3),CI(i,1)) = Count(CI(i,2),CI(i,3),CI(i,1)) + 1;
end

zDisplayPairCounts(Count,CountFilename);