% zBestRotation(X,Y) finds the least squares rotation
% of points X onto points Y
%
% X and Y are n by 3 matrices containing the locations of corresponding points
% What is returned is the best fit to Y = X*R'
% R is the 3x3 rotation matrix

% It is assumed that the means have been subtracted from X and Y

function [R] = zBestRotation(A,B)

n = length(A(:,1));                 % number of data points

M = zeros(3,3);                     % initialize the matrix M

for i=1:n,                          % go through all data points
  M = M + A(i,:)' * B(i,:);         % "outer product" of two vectors -> matrix
end

[v,d] = eig(M'*M);                  % eigenvector matrix, eigenvalue matrix
[e,i] = sort(-diag(d));             % sort the eigenvalues in decreasing order

a = v(:,i(1));                      % eigenvector with largest eigenvalue
b = v(:,i(2));
c = v(:,i(3));

Ma = M*a/sqrt(-e(1));
Mb = M*b/sqrt(-e(2));

if det([a b c]) > 0,
  g = cross(Ma,Mb);
else 
  g = -cross(Ma,Mb);
end

R = [a b c]*[Ma Mb g]';

