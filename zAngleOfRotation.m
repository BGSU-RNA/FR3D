% zAngleOfRotation calculates the angle of rotation from a rotation matrix R
% A good online reference is:
% http://www.mathworks.com/access/helpdesk/help/toolbox/physmod/mech/mech_review7.html

function [alpha] = zAngleOfRotation(R)

alpha = 2*acos(min(1,sqrt(trace(R)+1)/2));

% occasionally trace(R) is slightly greater than 3, making alpha complex
%[trace(R)-3] 
%[ sqrt(trace(R)+1)/2 alpha]
