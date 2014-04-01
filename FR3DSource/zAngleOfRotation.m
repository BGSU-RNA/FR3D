% zAngleOfRotation calculates the angle of rotation from a rotation matrix R

function [alpha] = zAngleOfRotation(R)

alpha = 2*acos(min(1,sqrt(trace(R)+1)/2));

% occasionally trace(R) is slightly greater than 3, making alpha complex
%[trace(R)-3] 
%[ sqrt(trace(R)+1)/2 alpha]

% There is an easier formula for angle of rotation, below.

if 0 > 1,
	theta = 31*pi/180;
	R1 = [[cos(theta) sin(theta) 0]; [-sin(theta) cos(theta) 0]; [0 0 1]];
	zAngleOfRotation(R1)*180/pi
	acos((trace(R1)-1)/2)*180/pi

	theta2 = 17*pi/180;
	R2 = [[cos(theta2) 0 sin(theta2)]; [0 1 0]; [-sin(theta2) 0 cos(theta2)]];
	zAngleOfRotation(R2)*180/pi
	acos((trace(R2)-1)/2)*180/pi

	zAngleOfRotation(R2*R1)*180/pi
	acos((trace(R2*R1)-1)/2)*180/pi

	zAngleOfRotation(R1'*R2)*180/pi
	acos((trace(R1'*R2)-1)/2)*180/pi
end
