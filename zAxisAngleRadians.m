% zAxisAngleBasic computes the axis and angle of rotation in an
% orthogonal matrix R

function [axis, angle] = zAxisAngleBasic(R)

[v,d] = eig(R);                % get eigenvectors and values of R
[y,i] = sort(imag(diag(d)));   % sort by imaginary part of eigenvalues

axis  = v(:,i(2));             % eigenvector with eigenvalue 1, imag part 0

[y,i] = sort(abs(axis));       % look for largest entries in axis

b       =  zeros(3,1);         % b will be perpendicular to axis
b(i(2)) =  axis(i(3));         % use the largest entries in axis
b(i(3)) = -axis(i(2));

angle = acos((b'*R*b) / (b'*b)) * sign(det([b R*b axis]));     
                               % angle of rotation, signed by rotation
                               % sense relative to axis of rotation
