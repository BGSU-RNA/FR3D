% zSubFile extracts nucleotides specified in i

function [F] = zSubFile(File,i)

F = File;
F.NT = File.NT(i);
F.NumNT = length(i);
F.Edge = File.Edge(i,i);
F.Coplanar = File.Coplanar(i,i);
F.Covalent = File.Covalent(i,i);
F.Crossing = File.Crossing(i,i);
if F.NumNT > 1,
  F.BasePhosphate = File.BasePhosphate(i,i);
else
  F.BasePhosphate = [];
end
if F.NumNT > 1,
  F.BaseRibose = File.BaseRibose(i,i);
else
  F.BaseRibose = [];
end
F.Backbone = File.Backbone(i,i);
if F.NumNT > 1,
  F.Range = File.Range(i,i);
else
  F.Range = [];
end
F.Flank = File.Flank(i,i);
F.Redundant = File.Redundant(i,i);
