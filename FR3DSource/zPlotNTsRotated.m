
function [void] = zPlotNTsRotated(File,Indices,VP,R,S)

for k=1:length(Indices),                 % Loop through all nucleotides
  zPlotOneNTRotated(File.NT(Indices(k)),VP,R,S);
end
