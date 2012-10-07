% zListPairData(Pair,h) prints pair analysis data
% if h == 0, no header is printed; if h == 1, a header is printed

function [void] = zListPairData(Pair,h)

if (h == 1),
  fprintf(' Sh(1)  Sh(2)  Sh(3) ');
  fprintf('Norm(3) PlAng    Ang    Gap   Comp   Hydrogen angles\n');
end

fprintf('%6.2f %6.2f ', Pair.Displ, Pair.Normal(3));
fprintf('%6.1f %6.1f ', Pair.PlaneAng, Pair.Ang);
fprintf('%6.2f %6.2f %6.2f %6.2f', Pair.Gap, Pair.Class);

if length(Pair.Hydrogen) > 0,
  for i=1:length(Pair.Hydrogen),
    fprintf('%6.1f ', Pair.Hydrogen(i).Angle);
  end
  fprintf('\n');
else
  fprintf('\n');
end
