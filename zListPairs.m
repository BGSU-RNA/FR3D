% zListPairs lists pair data in columns

function [void] = zListPairs(File,SP,ListMode,VP)

if ~isfield(VP,'ListItems')
  ListItems = [1 22 2 3 4 5 6 9 10 11 12 13 14 16 17 21];
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

M = [];                                     % accumulate data for summary

for k=1:length(SP)
  f = SP(k).Filenum;
  p = File(f).Pair(SP(k).PairIndex);                 % Current pair
  n1 = File(f).NT(p.Base1Index);
  n2 = File(f).NT(p.Base2Index);

  a = [];

 if ListMode == 2,
  if (mod(k,30) == 1),
    fprintf('%s\n',Header);
  end

  fprintf('%3d', k);                                 % Pair number
  for i = 1:length(ListItems),
    switch abs(ListItems(i)),
        case  1, fprintf('%12s',File(f).Filename);
        case  2, fprintf('%2s %4s',n1.Base, n1.Number);
        case  3, fprintf('%2s %4s',n2.Base, n2.Number);
        case  4, fprintf('%7.2f',p.Displ(1));
                 a = [a p.Displ(1)];
        case  5, fprintf('%7.2f',p.Displ(2));
                 a = [a p.Displ(2)];
        case  6, fprintf('%7.2f',p.Displ(3));
                 a = [a p.Displ(3)];
        case  7, fprintf('%7.2f',p.Normal(1));
                 a = [a p.Normal(1)];
        case  8, fprintf('%7.2f',p.Normal(2));
                 a = [a p.Normal(2)];
        case  9, fprintf('%7.2f',p.Normal(3));
                 a = [a p.Normal(3)];
        case 10, fprintf('%6.1f',p.Ang);
                 a = [a p.Ang];
        case 11, fprintf('%8.2f',SP(k).C1pC1p);
                 a = [a SP(k).C1pC1p];
        case 12, fprintf('%7.2f',p.Gap);
                 a = [a p.Gap];
        case 13, fprintf('%8.2f',p.MinDist);
                 a = [a p.MinDist];
        case 14, fprintf('%8.2f',p.Class);
        case 15, fprintf('%7.2f',SP(k).HandClass);
        case 16, fprintf('%7.2f',p.Classes(1));
        case 17, fprintf('%7.2f',p.Distances(1));
        case 18, if length(p.Hydrogen) > 0,
                   for i=1:length(p.Hydrogen),
                     fprintf('%6.1f ', p.Hydrogen(i).Angle);
                   end
                 else
                   fprintf('                 ');
                 end
        case 19, fprintf('%6.1f',p.StackingOverlap);
        case 21, fprintf('%7.2f',SP(k).PairDisc);
        case 22, if ~isempty(File(f).Info.Resolution),
                   fprintf('%6.2f',File(f).Info.Resolution);
                 else
                   fprintf('      ');
                 end
    end
  end
  fprintf('\n')

  M = [M; a];
end


 if mod(k,30) == 0,
   fprintf('Press return to continue, enter q to quit ');
   inp = input('','s');
   if ~isempty(inp),
     break
   end
 end

end



if length(M(:,1)) == 1,
  M = [M; M];
end

return

  fprintf('\n');
  fprintf('                                Sh(1)  Sh(2)  Sh(3) ');
  fprintf('Norm(3)   Ang    Gap MinDist\n');
  fprintf('Maximum                        ');
  fprintf('%6.2f ', max(M));
  fprintf('\n');
  fprintf('95th percentile                ');
  fprintf('%6.2f ', prctile(M,95));
  fprintf('\n');
  fprintf('Median                         ');
  fprintf('%6.2f ', median(M));
  fprintf('\n'); 
  fprintf(' 5th percentile                ');
  fprintf('%6.2f ', prctile(M,5));
  fprintf('\n');
  fprintf('Minimum                        ');
  fprintf('%6.2f ', min(M));
  fprintf('\n');


