% xPrecomputeForDiscrepancy(c) computes a few things that are needed for
% the full discrepancy calculation.  c is a structured variable with fields
% NT, LocWeight, NumNT

function [c] = xPrecomputeForDiscrepancy(c)

c.Centers = cat(1,c.NT.Center);

if c.NumNT > 2,
  c.WeightedCenter          = c.LocWeight * c.Centers / sum(c.LocWeight);           % a single vector, the normalized mean of the centers
  c.CenteredCenters         = c.Centers-ones(c.NumNT,1)*c.WeightedCenter;           % subtract the single center from each center
elseif c.NumNT == 2,
  c.R  = c.NT(2).Rot' * c.NT(1).Rot;
  c.T1 = (c.NT(2).Center - c.NT(1).Center)*c.NT(1).Rot;
  c.T2 = (c.NT(1).Center - c.NT(2).Center)*c.NT(2).Rot;
end
