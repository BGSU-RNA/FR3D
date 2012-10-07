% xDiscrepancy(Model,Cand) calculates the discrepancy between Model and 
% Cand, which is an array of NT's.
% As soon as the discrepancy exceeds Model.RelCutoff, the calculation stops.
% The current sum is returned as a negative discrepancy, shifted to
% record when the calculation stopped.

% Model and Cand are two list of indices or nucleotide numbers

function [Disc,R] = xDiscrepancy(File1,Model,File2,Cand,LocationWeight,AngleWeight)

L = length(Cand);

if nargin < 5,
  LocationWeight = ones(1,length(Model));
  AngleWeight    = ones(1,length(Model));
else
  LocationWeight = L * LocationWeight / sum(LocationWeight);
end

if length(Model) ~= length(Cand),           % motif sizes must be the same
  Disc = [];
else

% if File1 is a text string (filename), load the file

if strcmp(class(File1),'char'),
  File1name = File1;
  File1 = zGetNTData(File1name,0);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Model),'char'),
  Model = {Model};
end

if strcmp(class(Model),'cell'),
  Model = File1.NT(zIndexLookup(File1,Model));
else
  Model = File1.NT(Model);
end

% if File2 is a text string (filename), load the file and display

if strcmp(class(File2),'char'),
  File2name = File2;
  File2 = zGetNTData(File2name,0);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(Cand),'char'),
  Cand = {Cand};
end

if strcmp(class(Cand),'cell'),
  Cand = File2.NT(zIndexLookup(File2,Cand));
else
  Cand = File2.NT(Cand);
end

% ------------------------------ Calculate discrepancy

if (L == 2),                                % two-nucleotide motif

  ModelR  = Model(2).Rot' * Model(1).Rot;
  ModelT1 = (Model(2).Center - Model(1).Center)*Model(1).Rot;
  ModelT2 = (Model(1).Center - Model(2).Center)*Model(2).Rot;

  C    = Cand(1).Rot' * Cand(2).Rot;
  r1   = ModelR * C;
  ang1 = zAngleOfRotation(r1);
  r2   = ModelR' * C';
  ang2 = zAngleOfRotation(r2);

  t1   = ModelT1 - (Cand(2).Center - Cand(1).Center)*Cand(1).Rot;
  t2   = ModelT2 - (Cand(1).Center - Cand(2).Center)*Cand(2).Rot;

  v    = AngleWeight(1);

  Disc = (sqrt(t1*t1' + (v^2)*ang1^2) + sqrt(t2*t2' + (v^2)*ang2^2))/4;

  R = eye(2);                               % not really how to rotate them

else                                        % more than two nucleotides

  ModelCenters = cat(1,Model.Center);
  ModelWeightedCenter = LocationWeight * ModelCenters / L;
  MCC                 = ModelCenters-ones(L,1)*ModelWeightedCenter;

  CandiCenters = cat(1,Cand.Center);
  CMean = LocationWeight * CandiCenters / L;
  CC = CandiCenters - ones(L,1) * CMean;  % subtract mean
  
  R = zBestRotation(CC, diag(LocationWeight)*MCC);      % candidate onto model
  
  S = LocationWeight * sum(((MCC - CC*R').^2)')';  % distances between centers

  n = 1;                                    % nucleotide number for angles
  v = 4 * AngleWeight.^2;                   % precompute a little
  
  while (n <= L),
    angbytwo = acos(min(1,sqrt(trace(R*Cand(n).Rot*(Model(n).Rot)')+1)/2));
    S   = S + (angbytwo^2)*v(n);
    n   = n + 1;
  end

  Disc = sqrt(S)/L;

end
end