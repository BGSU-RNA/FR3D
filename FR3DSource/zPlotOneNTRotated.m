% zPlotOneNTRotated(NT,ViewParam,R,S) rotates NT and displays with ViewParam

function [void] = zPlotOneNTRotated(NT,ViewParam,R,S)

NewNT.Code  = NT.Code;
L           = length(NT.Fit(:,1));         % Number of base atoms

NewNT.Fit   = (NT.Fit   - ones(L,1) *S) * R; % rotated into position
s           = length(NT.Sugar(:,1));                % this sometimes varies
NewNT.Sugar = (NT.Sugar - ones(s,1)*S) * R; % rotate sugar too
NewNT.Base  = NT.Base;
NewNT.Number= NT.Number;

zPlotOneNT(NewNT,ViewParam)
