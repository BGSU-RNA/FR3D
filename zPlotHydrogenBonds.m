% zPlotHydrogenBonds(File,PairIndex,R,S) draws a dotted line between bases
% corresponding to a hydrogen bond

function [void] = zPlotHydrogenBonds(NT1,NT2,Class,R,S)

if (abs(Class) < 14),
  [X,Y] = zHydrogenLocations(NT1,NT2,Class);
else
  X = [];
  Y = [];
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


