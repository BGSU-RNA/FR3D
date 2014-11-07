% zPlotOneAARotated(AA,ViewParam,R,S) rotates AA and displays with ViewParam

function [void] = zPlotOneAARotated(AA,ViewParam,R,S)

NewAA = AA;

L           = length(AA.Loc(:,1));         % Number of AA atoms
NewAA.Loc   = (AA.Loc   - ones(L,1) *S) * R; % rotated into position

zPlotOneAA(NewAA,ViewParam)
