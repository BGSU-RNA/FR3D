% zPlotOneNTRotated(NT,ViewParam,R,S) rotates NT and displays with ViewParam

function [void] = zPlotOneNTRotated(NT,ViewParam,R,S)

NewNT.Code  = NT.Code;
L           = length(NT.Fit(:,1));         % Number of base atoms
NewNT.Fit   = (NT.Fit   - ones(L,1) *S) * R; % shift and rotate
NewNT.Base  = NT.Base;
NewNT.Number= NT.Number;

if NT.Code < 9                             % standard RNA or DNA
    s           = length(NT.Sugar(:,1));        % this sometimes varies
    NewNT.Sugar = (NT.Sugar - ones(s,1)*S) * R; % shift and rotate sugar
end

if NT.Code == 9,                  % observed atoms of modified nucleotide
    LocDict   = containers.Map;   % use a dictionary to store locations by atom name
    SugarDict = containers.Map;   % use a dictionary to store locations by atom name
    BetaDict  = containers.Map;   % use a dictionary to store locations by atom name
    for atom in NT.Loc
        LocDict(atom) = (NT.Loc(atom) - S) * R; % shift and rotate base atom
        SugarDict(atom) = (NT.Sugar(atom) - S) * R; shift and rotate sugar atom
    end
    NewNT.Loc = LocDict;
    NewNT.Sugar = SugarDict;
end

if isfield(NT,'Beta'),
  NewNT.Beta  = NT.Beta;
end

zPlotOneNT(NewNT,ViewParam)
