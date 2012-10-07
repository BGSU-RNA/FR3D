% zPairDiscrepancy(Pair1,Pair2) calculates a measure of the discrepancy
% between Pair1 and Pair2.

% This is a faulty calculation!
% Pair data alone doesn't give enough data to compute the FR3D discrepancy!

% What is below seems to be symmetric at least

function [d] = zPairDiscrepancy(Pair1,Pair2)

t1 = Pair1.Displ - Pair2.Displ;
r = Pair1.Rot * Pair2.Rot';

t2 = Pair1.Displ * Pair1.Rot - Pair2.Displ * Pair2.Rot;

[ax,ang] = zAxisAngleRadians(r);              % rotation angle without a flip

d = norm(t1) + norm(t2) + abs(ang)*(57.29577951308232/20); % convert to degrees, combine

%[norm(t) abs(ang)/20]
