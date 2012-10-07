% zBasePhosphateText(e) converts internal code e for a base-phosphate interaction to text for human use.  If edge = 1, it lists the active edge as well

function [EE,Base] = zBasePhosphateText(e,edge)

if nargin < 2,
  edge = 0;
end

E = [];
EE = [];
Base = [];

% The following vector converts internal FR3D codes to Base-phosphate categories0 to 9.  To change the categories, change here and in xGetEdgeNums.m

BPCat = [2 6 7 0 6 7 8 9 0 1 3 4 5 0 5 9 0 7 4];  % updated 8-19-2008
ae    = {'SW',' W',' H',' H',' W',' H',' H',' H',' H',' S',' W',' W',' W',' H',' W',' H',' H',' H',' W'};
Bases = {'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'G', 'G', 'U', 'U', 'U', 'C', 'G'};

for i=1:length(e),
  if e(i) > 100,
    E = [E 'n'];
    a = e(i) - 100;
  elseif e(i) < -100,
    E = [E 'n'];
    a = e(i) + 100;
  else
    E = [E ' '];
    a = e(i);
  end

  if abs(a) >= 1 && abs(a) <= 19,

    if a > 0,
      E = [E num2str(BPCat(a)) 'BPh'];
    elseif a < 0,
      E = [E num2str(BPCat(-a)) 'PhB'];
    else
      E = [E '  -'];
    end

    a = abs(fix(a));

    if edge > 0,
      if E(1) == ' ',
        E(1) = '_';
      end
      E = [ae{a} E];
    end

    EE = [EE E];                              % accumulate text
    Base = [Base Bases{a}];
  end
end
