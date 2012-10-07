% dNeedlemanWunsch(seq1,seq2,p,d) aligns sequences seq1 and seq2 using probability p of base conservation and gap penalty d.  align1 and align2 list which positions in seq1 and seq2 are aligned.  s1 and s2 are gapped versions of seq1 and seq2 which display the alignment.

function [matches,align1,align2,s1,s2] = dNeedlemanWunsch(seq1,seq2,p,d)

if nargin < 4,
  d = 2;                            % default gap penalty
end

if nargin < 3,
  p = 0.95;                         % default probability of base conservation
end

N = length(seq1);
M = length(seq2);

align1 = [];
align2 = [];

pa=1/N;
pb=1/M;

nwmatrix = zeros(N+1,M+1);          % indices are 1 to N+1
tracematrix = zeros(N+1,M+1);

for i = 1:N,
  nwmatrix(i+1,1) = -i*d;           % gap penalties along edge
  tracematrix(i+1,1) = 2;           % traceback option #2
end

for i = 1:M,
  nwmatrix(1,i+1) = -i*d;           % gap penalties along edge
  tracematrix(1,i+1) = 1;           % traceback option #1
end

S1 = seq1'*ones(1,M);               % spread characters into NxM matrix
S2 = ones(N,1)*seq2;                % spread characters into NxM matrix

S = log(p)*(S1==S2) + log((1-p))*(S1~=S2) - log(pa*pb);

for i = 1:N,                    % loop through all ways of extending alignment
  for j = 1:M,
    [a,b] = max([nwmatrix(i+1,j)-d nwmatrix(i,j+1)-d nwmatrix(i,j)+S(i,j)]);
    nwmatrix(i+1,j+1) = a;                % max value
    tracematrix(i+1,j+1) = b;             % which direction the max is from
  end
%  [a,b] = max(nwmatrix(i+1,:));          % doesn't do anything?
end

score = nwmatrix(N,M);

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

[m,a,b,s,t] = dNeedlemanWunsch(seq1,seq2,p,d)
[s;t]

seq1 = 'aacguuguggaa';
seq2 = 'aacguaugugcaga';
[Score, Alignment, Start] = swalign(seq1,seq2,'Alphabet','NT')

[Score, Alignment, Start] = swalign('uuuuaacguuguggaa','aacguaugugcaga','Alphabet','NT')

File = zAddNTData({'1s72','2aw4'});
seq1 = cat(2,File(1).NT.Base);
seq2 = cat(2,File(2).NT.Base);
[Score, Alignment, Start] = nwalign(seq1,seq2,'Alphabet','NT')
