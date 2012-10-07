% zPlotHydrogen(File,PairIndex,R,S) draws a dotted line between bases
% corresponding to a hydrogen bond

function [void] = zPlotHydrogen(File,PairIndex,R,S)

Pair  = File.Pair(PairIndex);
p     = Pair.Base1Index;
q     = Pair.Base2Index;
NT1   = File.NT(p);
NT2   = File.NT(q);
Class = Pair.Class;

if (abs(Class) < 14),
  [X,Y] = zHydrogenLocations(NT1,NT2,Class);
else
  ci = File.CI(p,q);                           % Comment index
  if ci > 0,                                   % Pair is in hand file
    HandClass = File.HandClass(ci);            % Use the hand class
  else
    HandClass = 0;                             % Hand class 0
  end
  [X,Y] = zHydrogenLocations(NT1,NT2,HandClass);
end

if ~isempty(X),
  L = length(X(:,1));
  A = (X - ones(L,1)*S) * R;                   % translate and rotate
  B = (Y - ones(L,1)*S) * R;                   % translate and rotate
  for i=1:L,
    C = [A(i,:); B(i,:)];
    plot3(C(:,1),C(:,2),C(:,3),':','Color','k','LineWidth',2);
  end
end


