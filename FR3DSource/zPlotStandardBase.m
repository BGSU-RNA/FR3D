% zPlotStandardBase(code,textoption) plots base with code,
% where A=1, C=2, G=3, U=4 and textoption:
% textoption = 0 - no text
% textoption = 1 - label atoms

% centeroption controls where the base is centered
% centeroption 0 - glycosidic atom at the origin (default)
% centeroption 1 - geometric center at the origin

function [void] = zPlotStandardBase(code,textoption,centeroption)

if nargin < 2,
  textoption = 0;
end

if nargin < 3,
  centeroption = 0;
end

  zStandardBases

  VP.Sugar = 0;

  BaseNames = 'ACGU';

  L = Lim(2,code);
  Q = StandardLoc(1:L,:,code);

  if centeroption == 1,
    M = Lim(1,code);
    Q = Q - ones(L,1)*mean(Q(1:M,:));
  end

  NT.Code = code;
  NT.Fit = Q;
  zPlotOneNT(NT,VP);

  H = [13 9 14 10];

  Z = [Q(H(code),:); Q(1,:)];
  k = [1 2]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',0.5*[1 1 1],'LineWidth',2,'LineStyle','-');
  if textoption > 0,
    hold on

    for j=1:L,
      text(Q(j,1), Q(j,2), 0.1, AtomNames{j,code},'FontSize',6);
    end
    title(['Standard ' BaseNames(code)])
  end

%  hold on
%  i = convhull(Q(:,1),Q(:,2));
%  plot(Q(i,1),Q(i,2),'k');

%  axis equal
