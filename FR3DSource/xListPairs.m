% zListPairs lists pair data in columns

function [void] = zListPairs(Pair,VP)

if ~isfield(VP,'ListItems')
  ListItems = [1 22 2 3 4 5 6 9 10 11 12 13 14 16 17 18];
else
  ListItems = VP.ListItems;
end

Header = ['    '];
for i = 1:length(ListItems),
  switch abs(ListItems(i)),
    case  1, Header = [Header '   Filename'];
    case  2, Header = [Header '  Nucl1'];
    case  3, Header = [Header '  Nucl2'];
    case  4, Header = [Header '  Disp1'];
    case  5, Header = [Header '  Disp2'];
    case  6, Header = [Header '  Disp3'];
    case  7, Header = [Header '  Norm1'];
    case  8, Header = [Header '  Norm2'];
    case  9, Header = [Header '  Norm3'];
    case 10, Header = [Header '   Ang'];
    case 11, Header = [Header ' C1*-C1*'];
    case 12, Header = [Header '    Gap'];
    case 13, Header = [Header ' MinDist'];
    case 14, Header = [Header '   Class'];
    case 15, Header = [Header '   Hand'];
    case 16, Header = [Header ' N.Exem'];
    case 17, Header = [Header '  EDist'];
    case 18, Header = [Header ' Hydrogen angles'];
    case 19, Header = [Header ' Overlap'];
    case 21, Header = [Header 'PairDisc'];
    case 22, Header = [Header ' Resol'];
  end
end

Text{1} = Header;

for k=1:length(Pair)
  p  = Pair(k);                 % Current pair
  n1 = p.NT1;
  n2 = p.NT2;

  Text{k+1} = '';

  Text{k+1} = [Text{k+1} sprintf('%3d', k)];                  % Pair number
  for i = 1:length(ListItems),
    switch abs(ListItems(i)),
        case  1, Text{k+1} = [Text{k+1} sprintf('%12s',p.Filename)];
        case  2, Text{k+1} = [Text{k+1} sprintf('%2s %4s',n1.Base, n1.Number)];
        case  3, Text{k+1} = [Text{k+1} sprintf('%2s %4s',n2.Base, n2.Number)];
        case  4, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Displ(1))];
        case  5, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Displ(2))];
        case  6, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Displ(3))];
        case  7, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Normal(1))];
        case  8, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Normal(2))];
        case  9, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Normal(3))];
        case 10, Text{k+1} = [Text{k+1} sprintf('%6.1f',p.Ang)];
        case 11, Text{k+1} = [Text{k+1} sprintf('%8.2f',p.C1pC1p)];
        case 12, Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Gap)];
        case 13, Text{k+1} = [Text{k+1} sprintf('%8.2f',p.MinDist)];
        case 14, Text{k+1} = [Text{k+1} sprintf('%8.2f',p.Class)];
        case 16, 
          if p.Distances(1) > 999,
            Text{k+1} = [Text{k+1} '  None'];
          else
            Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Classes(1))];
          end
        case 17, 
          if p.Distances(1) > 999,
            Text{k+1} = [Text{k+1} '   '];
          else
            Text{k+1} = [Text{k+1} sprintf('%7.2f',p.Distances(1))];
          end
        case 18, if length(p.Hydrogen) > 0,
                   for i=1:length(p.Hydrogen),
                     Text{k+1} = [Text{k+1} sprintf('%6.1f ', p.Hydrogen(i).Angle)];
                   end
                 else
                   Text{k+1} = [Text{k+1} sprintf('                 ')];
                 end
        case 19, Text{k+1} = [Text{k+1} sprintf('%6.1f',p.StackingOverlap)];
        case 22, Text{k+1} = [Text{k+1} sprintf('%6.2f',p.Resol)];
    end
  end
end

if isdeployed,
  mEditbox(Text,'List of pair parameters',10);
else
  for i=1:length(Text),
    fprintf('%s\n',Text{i});
  end
end



