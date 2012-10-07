% zAxisAngle(R) computes the axis and angle of rotation in an orthogonal
% matrix R, then makes the axis point up and the angle be between -90 and 270
% degrees


function [axis, angle] = zAxisAngle(R)

[axis, angle] = zAxisAngleRadians(R);

angle = angle*57.29577951308232;

if axis(3) < 0,                   % make axis point up
  axis = -axis;
  angle = -angle;
end

if angle < -90,	          % make angles be between -90 and 270 degrees
  angle = angle + 360;
end

