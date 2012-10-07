% zFileNumberLookup(File,Name) finds the index in the array File of the molecule called Name


function [f] = zFileNumberLookup(File,Name)

f = 0;

for i=1:length(File),
  if strcmp(File(i).Filename,Name),
    f = i;
  end
end
