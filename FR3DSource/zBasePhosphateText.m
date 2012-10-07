
function [E] = zBasePhosphateText(e)

E = [];

% The following vector converts internal FR3D codes to Base-phosphate categories1 to 9.  To change the categories, change here and in xGetEdgeNums.m

BPCat = [2 5 6 9 5 6 7 8 9 1 3 4 5 9 4 8 9];

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

  if a > 0,
    E = [E num2str(BPCat(a)) 'BP'];
  elseif a < 0,
    E = [E num2str(BPCat(-a)) 'PB'];
  else
    E = [E '  -'];
  end
end
