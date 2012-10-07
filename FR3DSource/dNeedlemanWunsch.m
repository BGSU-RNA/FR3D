


function [matches,align1,align2,s1,s2] = NeedlemanWunsch(seq1,seq2,p,d)

N = length(seq1);
M = length(seq2);

align1 = [];
align2 = [];

pa=1/N;
pb=1/M;

nwmatrix = zeros(N+1,M+1);          % indices are 1 to N+1
tracematrix = zeros(N+1,M+1);

for i = 1:N,
  nwmatrix(i+1,1) = -i*d;
  tracematrix(i+1,1) = 2;
end

for i = 1:M,
  nwmatrix(1,i+1) = -i*d;
  tracematrix(1,i+1) = 1;
end

for i = 1:N,
  for j = 1:M,
    if seq1(i) == seq2(j),
      pab = p;
    else
      pab = 1-p;
    end
    s = log(pab) - log(pa*pb);
    [a,b] = max([nwmatrix(i+1,j)-d nwmatrix(i,j+1)-d nwmatrix(i,j)+s]);
    nwmatrix(i+1,j+1) = a;                % max value
    tracematrix(i+1,j+1) = b;             % which direction the max is from
  end
  [a,b] = max(nwmatrix(i+1,:));
end

%round(nwmatrix)
%tracematrix

% output

i = N;     % current location in seq1
j = M;     % current location in seq2
s1 = [];   % accumulating letters in first row of alignment
s2 = [];   % accumulating letters in second row of alignment

matches = 0;

while (i > 0) || (j > 0),
  switch tracematrix(i+1,j+1),
  case 1,
    s1 = ['-' s1];
    s2 = [seq2(j) s2];
    j = j - 1;
  case 2,
    s1 = [seq1(i) s1];
    s2 = ['-' s2];
    i = i - 1;
  case 3,
    s1 = [seq1(i) s1];
    s2 = [seq2(j) s2];
    align1 = [i align1];
    align2 = [j align2];
    i = i - 1;
    j = j - 1;
    matches = matches + 1;
  end
end    



return

seq1 = 'avseryu';
seq2 = 'atvpeymuq';

seq1 = 'abczzdefghi';
seq2 = 'wwwwabcdefghijklm';

p = .995;                        % probability an A stays an A
d = 2;                         % gap penalty

[m,a,b,s,t] = nw(seq1,seq2,p,d)
[s;t]
