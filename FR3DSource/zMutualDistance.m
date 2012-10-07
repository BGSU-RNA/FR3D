% zMutalDistance(A,L) finds the mutual distances between the rows of A and
% returns a sparse matrix D in which all entries are less than L

function [D] = zMutualDistance(A,L)

[s,t] = size(A);                       % s is the number of rows

if s < 2000,                           % calculate D all at once
  D = zDistance(A);
  D = sparse(D .* (D < L));

else                                   % calculate small parts of D at a time
  if s < 6000,                         % decide on the size of these blocks
    N = ceil(s/3);                     % s/3 is optimal
  elseif s < 8000,
    N = ceil(s/4);                     % s/4 is also good, keeps N small
  elseif s < 10000,
    N = ceil(s/5);
  else
    N = 2000;                          % no point going above 2000
  end

  u = 1;                               % lower limit of the block
  w = N;                               % upper limit
  c = 1;                               % block counter
  g = 2;                               % keep going until g = 0
  while g > 0,
    i{c} = (u:w)';                     % form a new block
    c = c + 1;
    u = u + N;                         % move to the next block
    w = w + N;
    if w > s && g == 2,                % first time above s
      w = s;
      g = 1;
    elseif w > s,
      w = s;
      g = 0;
    end
  end
  
  D = sparse(s,s);                     % sparse matrix to store dists

  mm = [];
  nn = [];
  vv = [];

  for a = 1:length(i),
    for b = a:length(i),

      DD = zDistance(A(i{a},:),A(i{b},:));
      DD = sparse(DD .* (DD < L));

      [m,n,v] = find(DD);              % non-zero entries

      if ~isempty(m),
        mm = [mm; i{a}(m)];
        nn = [nn; i{b}(n)];
        vv = [vv; v];

        if b > a,
          mm = [mm; i{b}(n)];            % symmetrize
          nn = [nn; i{a}(m)];
          vv = [vv; v];
        end
      end

    end
  end

  D = sparse(mm,nn,vv,s,s);
end
