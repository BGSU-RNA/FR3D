% zPlotOneAARotated(AA,ViewParam,R,S) rotates AA and displays with ViewParam

function [void] = zPlotOneAARotated(AA,ViewParam,R,S)

L           = length(AA.Loc(:,1));         % Number of AA atoms
NewAA.Loc   = (AA.Loc   - ones(L,1) *S) * R; % rotated into position
NewAA.Unit  = AA.Unit;
NewAA.Number= AA.Number;

if isfield(AA,'Beta'),
  NewAA.Beta  = AA.Beta;
end

zPlotOneAA(NewAA,ViewParam)
