% xSuperimposeCandidates(Model,Cand) returns the rotation matrix R and
% shift vector S that will align the centers of nucleotides in Cand with
% those in Model.  It also returns E, the least-squares error.

function [R,S,E] = xSuperimposeCandidates(Model,Cand,LocationWeight,AngleWeight)

L = length(Model);

if nargin < 5,
  LocationWeight = ones(1,length(Model));
  AngleWeight    = ones(1,length(Model));
else
  LocationWeight = L * LocationWeight / sum(LocationWeight);
end

if length(Model) ~= length(Cand),           % motif sizes must be the same
  R = [];
  S = [];
  E = [];
else

% ------------------------------

  ModelCenters        = cat(1,Model.Center);
  ModelWeightedCenter = LocationWeight * ModelCenters / L;
  MCC                 = ModelCenters-ones(L,1)*ModelWeightedCenter;

  CandCenters         = cat(1,Cand.Center);
  CandWeightedCenter  = LocationWeight * CandCenters / L;
  CC                  = CandCenters - ones(L,1) * CandWeightedCenter;

  if (L == 2),                                % two-nucleotide motif

    R = Model(1).Rot * Cand(1).Rot';
    S = CandWeightedCenter;
    E = [];

  else                                        % more than two nucleotides
  
    R = zBestRotation(CC, diag(LocationWeight)*MCC);% candidate onto model
    S = CandWeightedCenter;
    E = LocationWeight * sum(((MCC - CC*R').^2)')'; % distances between centers

  end
end