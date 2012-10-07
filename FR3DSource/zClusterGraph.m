% zClusterGraph(D,Lab,W) treats D as the distances between instances with labels given by Lab.  It re-orders the instances to group them into clusters, and puts nearby clusters next to one another.  Then it displays the distances graphically.  W is the number of characters of Lab to use.  pp is an optional user-supplied ordering vector.

function [p] = zClusterGraph(D,Lab,W,pp,Table)

[s,t] = size(D);

NoLabel = 0;

if (nargin < 2) || (length(Lab) == 0),
  for i = 1:s,
    Lab{i} = '    ';
  end
  W = 1;
  NoLabel = 1;
end  

if nargin < 3,
  W = length(Lab{1});
  for i = 1:s,
    W = min(length(Lab{i}),W);
  end
end  

if nargin < 4,
  pp = [];
end

if nargin < 5,
  Table = 0;
end

if isempty(pp),
  p = zOrderbySimilarity(D);              % put instances in a useful order
else
  p = pp;                                 % use user-supplied ordering
end

DDD = D(p,p);                            % re-order according to p

% ------------------------------------------ Print table

if Table == 1,
  fprintf('                             ');
  for j=1:t,
    w = min(length(Lab{p(j)}), 6);
    fprintf('%6s', Lab{p(j)}(1:w));
  end
  fprintf('\n');
  for i=1:s,
    fprintf('%24s ', Lab{p(i)}(1:28));
    for j= 1:t,
      fprintf('%5.2f ', DDD(i,j));
    end
    fprintf('\n');
  end
  fprintf('\n');
end

% ------------------------------------------ Display graph of isodiscrepancies

DDDD = zeros(s+1,t+1);           % pad with zeros
DDDD(1:s,1:t) = DDD;

pcolor(real(DDDD))

%pcolor(2*floor(DDDD/2))

shading flat
axis ij
view(2)
%colormap('default');
%map = colormap;
%map = map((end-8):-1:8,:);
%map = map((end-8):-1:end,:);
%colormap(map);
%caxis([0 16]);
%colorbar('location','eastoutside');

if Table > 0,
  for i = 1:s,
    for j = 1:s,
      tt = sprintf('%3.1f',DDDD(i,j));
      text(i+0.5,j+0.5,tt,'HorizontalAlignment','Center');
    end
  end
end

for i = 1:t,
  SLab{i} = Lab{i}(1:W(1));
end

if length(W) > 1,
  for i = 1:t,
    SSLab{i} = Lab{i}(1:W(2));
  end
end

if (length(Lab) < 20) || (NoLabel == 1),
  FS = 12;
elseif length(Lab) < 50,
  FS = 8;
elseif length(Lab) < 80,
  FS = 4;
elseif length(Lab) < 100,
  FS = 3;
elseif length(Lab) < 300,
  FS = 2;
elseif length(Lab) < 400,
  FS = 1;
else
  FS = 0.5;
end

if length(W) > 1,
  set(gca,'XTick',(1:s)+0.5)
  set(gca,'XTickLabel',SSLab(p),'FontSize',FS)
else
  XFS = min(6,FS);
  for i = 1:t,
    text(i+0.5,t+1,SLab{p(i)},'Rotation',90,'HorizontalAlignment','right','FontSize',XFS);
  end
  set(gca,'XTickLabel',[])
end

set(gca,'YTick',(1:s)+0.5)
set(gca,'YTickLabel',SLab(p),'FontSize',FS)

%set(gca,'OuterPosition',[0.1 0.1 1.1 1.1]); % doesn't help with pdf file!
