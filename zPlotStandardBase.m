% zPlotStandardBase(code,textoption) plots base with code,
% where A=1, C=2, G=3, U=4 and textoption:
% textoption = 0 - no text
% textoption = 1 - label atoms

function [void] = zPlotStandardBase(code,textoption)

if nargin < 2,
  textoption = 0;
end

  zStandardBases

  VP.Sugar = 0;

  BaseNames = 'ACGU';

  if ~exist('Centered'),
    Centered = 0;
  end

  L = Lim(2,code);
  Q = StandardLoc(1:L,:,code);

  if Centered == 1,
    M = Lim(1,code);
    Q = Q - ones(L,1)*mean(Q(1:M,:));
  end

  NT.Code = code;
  NT.Fit = Q;
  zPlotOneNT(NT,VP);

  if textoption > 0,
    hold on

    for j=1:L,
      text(Q(j,1)+0.1, Q(j,2), AtomNames{j,code});
    end
    title(['Standard ' BaseNames(code)])
  end

%  hold on
%  i = convhull(Q(:,1),Q(:,2));
%  plot(Q(i,1),Q(i,2),'k');

%  axis equal
