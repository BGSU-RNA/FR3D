% zNeedlemanWunschAffineGap aligns sequences S and T using the GapOpen and GapExtend scores, which should be specified as positive numbers.
% If the letters in S and T are ACGU or ACGT, Matlab's nuc44 scoring is used for substitutions; 5 for match, -4 for any substitution.
% If the letters differ from this, the match score is 1 and the substitution score is -Inf, so substitutions are not allowed.

function [matches,align1,align2,s1,s2] = zNeedlemanWunschAffineGap(S,T,GapOpen,GapExtend,MatchScore,LetterKey)

if nargin < 1,
  S = 'ACGUACGUACGU';
  T = 'ACUACCGCCCAGU';
  T = 'ACGUCCGUACGU';
  T = 'ACGUUGCAUGCAACGUACGU';
end

if nargin < 5,
  LetterKey = unique([S T]);
  L = length(LetterKey);
  if strcmp(LetterKey,'ACGU') || strcmp(LetterKey,'ACGT'),
    match = 5;
    ti = -4;     % transition
    tv = -4;     % transversion
    MatchScore = [match tv ti tv; tv match tv ti; ti tv match tv; tv ti tv match];
  else
    MatchScore = -Inf * ones(L,L);           % only allow exact matches
    for i = 1:L,
      MatchScore(i,i) = 1;
    end
  end
end

if nargin < 3,
  GapOpen   = -8;
  GapExtend = -4;
else
  GapOpen = -GapOpen;
  GapExtend = -GapExtend;
end

% ------------------------- convert sequences to numeric

clear SS TT

for i = 1:L,
  j = find(S == LetterKey(i));
  SS(j) = i;
  j = find(T == LetterKey(i));
  TT(j) = i;
end

A = length(SS);
B = length(TT);

% -------------------------- Fill in dynamic programming matrices and pointers

Max        = zeros(A+1,B+1);
MaxAbove   = zeros(A+1,B+1);
MaxLeft    = zeros(A+1,B+1);
Point      = zeros(A+1,B+1);

% Point is the pointer for Max
%   1 means the max came from extending a gap from above
%   2 means the max came from opening a gap from above
%   3 means the max came from matching two letters
%   4 means the max came from opening a gap from the left
%   5 means the max came from extending a gap from the left

% ----------------- Fill in aligning S to gaps
for i = 2:(A+1),
  Max(i,1)        = GapOpen + (i-2) * GapExtend;
  Point(i,1)      = 1;
  MaxAbove(i,1)   = GapOpen + (i-2) * GapExtend;
  MaxLeft(i,1)    = -Inf;
end

Point(2,1) = 2;

% ----------------- Fill in aligning T to gaps
for j = 2:(B+1),
  Max(1,j)        = GapOpen + (j-2) * GapExtend;
  Point(1,j)      = 5;
  MaxLeft(1,j)    = GapOpen + (j-2) * GapExtend;
  MaxAbove(1,j)   = -Inf;
end

Point(1,2) = 4;

% ----------------- Fill in the rest of the matrices using recursions
for i = 2:(A+1),
  for j = 2:(B+1),
    [y,p] = max([(MaxAbove(i-1,j) + GapExtend) (Max(i-1,j) + GapOpen)]);
    MaxAbove(i,j)   = y;

    [y,p] = max([(MaxLeft(i,j-1) + GapExtend) (Max(i,j-1) + GapOpen)]);
    MaxLeft(i,j)    = y;

    [y,p] = max([(MaxAbove(i-1,j) + GapExtend) (Max(i-1,j) + GapOpen) (Max(i-1,j-1)+MatchScore(SS(i-1),TT(j-1))) (Max(i,j-1) + GapOpen) (MaxLeft(i,j-1) + GapExtend)]);
    Max(i,j)        = y;
    Point(i,j)      = p;
  end
end

score = Max(A+1,B+1);

% ----------------- Traceback

i = A+1;
j = B+1;

align1 = [];
align2 = [];
s1     = '';
s2     = '';

while i > 1 || j > 1,

%  fprintf('Letter %3d of S, letter %3d of T, p = %d, i = %3d, j = %3d.\n', i-1, j-1, Point(i,j), i, j);

  switch Point(i,j),
  case 1


%Point(1:i,j)

%pause

    while Point(i,j) == 1,
      i  = i - 1;
      s1 = [S(i) s1];
      s2 = ['-' s2];
    end

    i  = i - 1;
    s1 = [S(i) s1];
    s2 = ['-' s2];

%[1 Point(i,j)]

%        fprintf('  MaxAbove pointer is p = %d\n', p);

  case 2
    i  = i - 1;
    s1 = [S(i) s1];
    s2 = ['-' s2];

%        fprintf('  MaxAbove pointer is p = %d\n', p);

  case 3                                         % max came from Max

%    fprintf('  Max pointer is p = %d\n', p);

    i      = i - 1;
    j      = j - 1;
    align1 = [i align1];
    align2 = [j align2];
    s1     = [S(i) s1];
    s2     = [T(j) s2];

  case 4
    j  = j - 1;
    s1 = ['-' s1];
    s2 = [T(j) s2];

  case 5
    while Point(i,j) == 5,
      j  = j - 1;
      s1 = ['-' s1];
      s2 = [T(j) s2];
    end

    j  = j - 1;
    s1 = ['-' s1];
    s2 = [T(j) s2];

  end
end

matches = sum(s1 == s2);
