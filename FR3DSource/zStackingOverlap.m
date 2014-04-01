% zStackingOverlap(N1,N2) computes a measure of overlap between nucleotides
% N1 and N2, by projecting N2 onto the plane of N1.

function [SO] = zStackingOverlap(N1,N2)

switch N1.Code,
  case 1,
X = [[ -2.100463  0.000000  4.271447  1.920945  0.230436 -2.100463]; ...
     [  0.447145 -1.009320  1.317924  5.150733  4.699718  0.447145]]';
t = [[    5    5    5]; ...
     [    2    4    2]; ...
     [    6    3    3]]';
  case 2,
X = [[ -2.082733  0.000000  2.269450  1.203833 -0.527970 -2.036772 -2.082733]; ...
     [  0.123632 -1.010259 -0.120783  4.411996  4.602202  2.647095  0.123632]]';
t = [[    6    6    4    4]; ...
     [    2    2    6    6]; ...
     [    7    3    3    5]]';
  case 3,
X = [[ -2.101572  0.000000  4.872516  5.295175  3.613335  1.396986 -0.751391 -2.101572]; ...
     [  0.463584 -1.009529  0.097781  1.782283  3.374547  4.394213  2.120132  0.463584]]';
t = [[    7    5    5    5    5]; ...
     [    2    7    3    7    3]; ...
     [    8    2    2    6    4]]';
  case 4,
X = [[ -2.082780  0.000000  2.292490  2.092152  0.177156 -2.124577 -2.082780]; ...
     [  0.111836 -1.008947 -0.048394  2.445179  4.020060  2.616537  0.111836]]';
t = [[    4    4    7    7]; ...
     [    5    2    4    4]; ...
     [    6    3    6    2]]';
end

S = N1.Fit(1,:);
R = N1.Rot;
L = length(N2.Fit(:,1));         % Number of base atoms
Y = (N2.Fit   - ones(L,1)*S) * R; % rotated into position

warning off

if exist('DelaunayTri') == 2,
  Y = Y(:,1:2);
  t = DelaunayTri(X);
  d = pointLocation(t,Y);
else
  d = tsearch(X(:,1),X(:,2),t,Y(:,1),Y(:,2)); % old method, no longer in Matlab
end

warning on

SO = length(find(d < 100));

%figure(2)
%clf
%plot(X(:,1),X(:,2),'k');
%hold on
%plot(Y(:,1),Y(:,2),'*');
