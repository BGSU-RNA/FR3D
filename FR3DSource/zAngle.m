% zAngle(A,B,C) computes the angle at B made by vectors BA and BC

function [theta] = zAngle(A,B,C)

d = A-B;
e = C-B;

theta = 180*acos(d*e'/(norm(d)*norm(e)))/pi;

