% zColorPairs(File,SP,Param,ViewParam) assigns color codes to each pair in
% File for use by zScatterPairs, according to the viewing parameters in
% ViewParam.


function [SP] = zColorPairs(File,SP,Param,ViewParam)

for k = 1:length(SP),                                     % Loop through pairs
  p  = File(SP(k).Filenum).Pair(SP(k).PairIndex);         % Current pair

  % -------- Round categories to nearest integer, if desired ------------- 

  if Param.Decimal == 1,                  % select with decimal places in categ
    pClass = p.Class;
    hClass = SP(k).HandClass; 
    eClass = p.Classes(1);
  else
    pClass = fix(p.Class);
    hClass = fix(SP(k).HandClass);
    eClass = fix(p.Classes(1));
  end

  switch ViewParam.Color,
    case 1, if any(pClass == Param.Category),
              SP(k).Color = 2;
            else
              SP(k).Color = 0;
            end
            if any(hClass == Param.Category),
              SP(k).Color = SP(k).Color + 4;
            end
    case {2, 3, 4}, SP(k).Color = p.Class;
    case 5, SP(k).Color = p.Paircode;
    case 6, SP(k).Color = p.Ang;
    case 7, SP(k).Color = p.Gap;
    case 8, if any(pClass == Param.Category),
              SP(k).Color = 2;
            else
              SP(k).Color = 0;
            end
            if any(eClass == Param.Category),
              SP(k).Color = SP(k).Color + 4;
            end
    case 9, SP(k).Color = p.Distances(1);
    case 10, SP(k).Color = p.Classes(1);
    case 11, SP(k).Color = p.StackingOverlap;
    case 12, SP(k).Color = p.Class;
    case 13, SP(k).Color = SP(k).PairDisc;
  end

end

