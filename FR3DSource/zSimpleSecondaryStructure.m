
function [C] = zSimpleSecondaryStructure(File)

N = File.NumNT;

md = 10;                                        % masking depth

A = fix(abs(File.Edge(1:N,1:N)));
C = A .* (A < 20);
C = triu(C);
M = triu(ones(size(C)));                       % full matrix

Orig = nonzeros(C);                            % original pair categories
fprintf('Originally %d basepairs\n', nnz(C));

figure(1)
clf
spy(C)
title(['Original basepairs in ' File.Filename]);

for p = 1:(N-1),                               % distance from diagonal

  k = find(diag(C,p)==1);                      % nonzero rows above diagonal

  for m = 1:length(k),
    i = k(m);                                  % row
    j = i + p;                                 % column

    if C(i,j) == 1,                            % cWW basepair between i and j
      d = md;
    else
      d = 2*md;                              % larger depth for non WC
    end

    a = 1:i;
    b = i:j;

    S = C(a,b);                              % store for later
    T = M(a,b);

    for x = 1:length(a),
      for y = 1:length(b),
        C(a(x),b(y)) = 0;
        M(a(x),b(y)) = 0;
      end
    end

    for x = 1:min(d,length(a)),
      for y = 1:min(d+1-x,length(b)),
        C(i-x+1,j-y+1) = S(end+1-x,end+1-y);
        M(i-x+1,j-y+1) = T(end+1-x,end+1-y);
      end
    end

    a = i:j;
    b = j:N;

    S = C(a,b);                              % store for later
    T = M(a,b);

    C(a,b) = zeros(size(C(a,b)));
    M(a,b) = zeros(size(M(a,b)));

    for x = 1:min(d,length(a)),
      for y = 1:min(d+1-x,length(b)),
        C(i+x-1,j+y-1) = S(x,y);
        M(i+x-1,j+y-1) = T(x,y);
      end
    end

  end

  fprintf('Diagonal %4d\n',p);
  if length(k) > 10000000,
    figure(1)
    clf
    spy(C);
    title(['Reduced basepairs in ' File.Filename]);
    figure(2)
    clf
    spy(M);
    title('Basepair mask to eliminate long-range non-nested basepairs');
    drawnow
  end
end

fprintf('At the end, %d basepairs\n', nnz(C));

Final = nonzeros(C);

for i = 1:15,
  fprintf('%4s original %3d final %3d\n', zEdgeText(i), length(find(Orig==i)), length(find(Final==i)));
end

C = triu(C) + triu(C)';

% Statistics for 2avy with depth 10:

%Originally  678 basepairs
%At the end, 577 basepairs
%cWw  original 472 final 460
%tWw  original   7 final   6
%cWH  original  12 final   5
%tWH  original  33 final  32
%cWS  original  12 final   8
%tWS  original   9 final   6
%cHH  original   0 final   0
%tHH  original   3 final   2
%cHS  original   7 final   4
%tHS  original  44 final  40
%cSs  original  45 final   7
%tSs  original  27 final   5
%bif  original   3 final   2
%Rib  original   4 final   0
%  -  original   0 final   0