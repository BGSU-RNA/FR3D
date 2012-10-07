% xDiscrepancy(Model,Cand) calculates the discrepancy between Model and 
% Cand, which is an array of NT's.
% One must take the square root and divide by the number of nucleotides after.
% As soon as the discrepancy exceeds Model.RelCutoff, the calculation stops.
% The current sum is returned as a negative discrepancy, shifted to
% record when the calculation stopped.

function [Disc] = xDiscrepancyFast(Model,Cand)

L = length(Cand);

if (L == 2),

  C    = Cand(1).Rot' * Cand(2).Rot;
  r1   = Model.R * C;
  ang1 = zAngleOfRotation(r1);
  r2   = Model.R' * C';
  ang2 = zAngleOfRotation(r2);

  t1   = Model.T1 - (Cand(2).Center - Cand(1).Center)*Cand(1).Rot;
  t2   = Model.T2 - (Cand(1).Center - Cand(2).Center)*Cand(2).Rot;

  v    = Model.AngleWeight(1);

  Disc = (sqrt(t1*t1' + (v^2)*ang1^2) + sqrt(t2*t2' + (v^2)*ang2^2))/4;

  if (Disc > Model.LDiscCutoff),
    Disc = -1;
  end

else

  MCC = Model.CenteredCenters;          % nucleotide centers in model
  MWC = Model.WeightedCenteredCenters;

  CandiCenters = cat(1,Cand.Center);
  CMean = Model.LocWeight * CandiCenters / Model.NumNT;

  CC = CandiCenters - ones(L,1) * CMean;  % subtract mean
  
  R = zBestRotation(CC, MWC);             % candidate onto model
  
  S = Model.LocWeight * sum(((MCC - CC*R').^2)')';  % distances between centers
  n = 1;                                    % nucleotide number for angles
  
  while (S <= Model.LDiscCutoff) & (n <= L),
    ang = zAngleOfRotation(R*Cand(n).Rot*(Model.NT(n).Rot)');
    S   = S + (ang^2)*Model.AngleWeight(n)^2;
    n   = n + 1;
  end

  if (n > 1) & (S <= Model.LDiscCutoff),
    Disc = S;
  else
    Disc = -1;                                % signal large discrepancy
  end

end

Disc = real(Disc);                            % sometimes not real, even 
