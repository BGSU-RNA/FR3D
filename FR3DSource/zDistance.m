% zDistance(A,B) finds the Euclidean distances between the rows of A and of B

function [D] = zDistance(A,B)

[M,t] = size(A);

if nargin == 2,
  [N,s] = size(B);
else
  N = 0;
  s = 0;
end

if nargin < 2,

  G = A * A';                        % inner products of rows of A

  if M < 150,
    X = diag(G) * ones(1,M);           % M by M matrix
  else
    X = repmat(diag(G),1,M);           % repeat a in each column
  end

  D = sqrt(X + X' - 2*G);            % |u-v| = sqrt(|u|^2 + |v|^2 - 2 u . v)

elseif s == t,

  a = sum((A.*A)');                   % sum of squares of each row of A
  b = sum((B.*B)');                   % sum of squares of each row of B

  if M < 150,
    X = a' * ones(1,N);                 % M by N matrix
    Y = ones(M,1) * b;                  % M by N matrix
  else
    X = repmat(b,M,1);                  % repeat b in each row
    Y = repmat(a',1,N);                 % repeat a in each column
  end

  D = sqrt(X + Y - 2*A*B');           % |u-v| = sqrt(|u|^2 + |v|^2 - 2 u . v)

else

  D = [];
  fprintf('zDistance: Matrix sizes are not compatible\n');

end
