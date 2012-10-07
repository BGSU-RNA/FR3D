% zSortPairs sorts according to the desired sort keys

function [SP] = zSortPairs(File,SP,ViewParam)

if length(ViewParam.SortKeys) > 0,

  M = [];                                     % accumulate data for sorting

  for k=1:length(SP)
    f = SP(k).Filenum;
    p = File(f).Pair(SP(k).PairIndex);                 % Current pair
    theta = -180*atan2(p.Displ(1),p.Displ(2))/pi;    
    if length(p.Hydrogen(:)) > 0,
      hydangles = cat(1,p.Hydrogen(:).Angle);
      hyddist   = cat(1,p.Hydrogen(:).Distance);
      b = [min(hydangles) max(hydangles) min(hyddist) max(hyddist)];
    else
      b = [200 0 0 0];
    end
    a = [];
    for i = 1:length(ViewParam.SortKeys),
      switch abs(ViewParam.SortKeys(i)),
        case  1, a = [a f];
        case  2, a = [a p.Paircode];
        case  3, a = [a p.Displ(1)];
        case  4, a = [a p.Displ(2)];
        case  5, a = [a p.Displ(3)];
        case  6, a = [a p.Normal(3)];
        case  7, a = [a p.PlaneAng];
        case  8, a = [a p.Ang];
        case  9, a = [a abs(p.Gap)];
        case 10, a = [a p.Class];
        case 11, a = [a SP(k).HandClass];
        case 12, a = [a SP(k).MinDist];
        case 13, a = [a p.Base1Index];
        case 14, a = [a min(p.Base1Index,p.Base2Index)];
        case 15, a = [a theta];
        case 16, a = [a b(1)];
        case 17, a = [a b(2)];
        case 18, a = [a b(3)];
        case 19, a = [a b(4)];
        case 20, a = [a p.Distances(1)];
        case 21, a = [a p.Classes(1)];
        case 22, a = [a p.StackingOverlap];
        case 23, a = [a SP(k).PairDisc];
        case 24, a = [a (abs(p.Gap)+abs(1-p.Normal(3))+SP(k).MinDist)];
      end
      if ViewParam.SortKeys(i) < 0,   % reverse order of sort
        a(i) = -a(i);
      end
    end
    M = [M; a];
  end

  [a,b] = sortrows(M);

  SP = SP(b);  
end