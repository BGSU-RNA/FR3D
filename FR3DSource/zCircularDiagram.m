% zCircularDiagram(File,Thickness) plots the pairwise interactions in File
% using colored chords around a unit circle.  Thickness controls the
% thickness of the lines, for different graphic output formats.

% zCircularDiagram('1s72') will load the file and display the diagram.

% Helpful suggestions for how to save the figure as a graphic file:
%  clf
%  zCircularDiagram(File(f),1);
%  saveas(gcf,[mypath FN '_circular_diagram.png'],'png');
%  [X,map] = imread([mypath FN '_circular_diagram.png']);
%  Y = X(30:830,210:1030,:);
%  imwrite(Y,[mypath FN '_circular_diagram.png']);

%  clf
%  zCircularDiagram(File(f),0.1);
%  saveas(gcf,[mypath FN '_circular_diagram.pdf'],'pdf');

function [void] = zCircularDiagram(File,Thickness)

if nargin < 2,
  Thickness = 1;
end

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

E  = fix(abs(File.Edge));
B  = E .* (E > 0) .* (E < 24);                 % pairs and stacks
R  = File.Range;

Edge = triu(B);                         % which pairs of bases to display

Color = (B==1).*(R==0) + 2*(B>1).*(B<14).*(R==0) + 3*(B==1).*(R>0) + 4*(B > 1).*(B < 14) .*(R>0) + 5*(B > 20) .* (B < 25);
                                        % disjoint possibilities

BP = abs(File.BasePhosphate);           % 
BP = (BP > 0) .* (BP < 100);            % exclude near BP and self interactions

[i,j,c] = find(triu(Color));
[ii,jj,cc] = find(6*BP);                % handle this separately, don't add
k = find(ii ~= jj);                     % eliminate self interactions

i = [i; ii(k)];
j = [j; jj(k)];
c = [c; cc(k)];

if length(i) > 0,

[s,t] = size(Edge);

% ---------------------------------- Draw the outside circle
% ---------------------------------- Later, break the circle for chains

tp = -6.2;                           % this value of 2pi leaves a gap
theta = 0:-0.01:tp;
plot(cos(theta), sin(theta),'k');
hold on

if s > 1000,
  step  = 50;
  sstep = 10;
elseif s > 500,
  step = 20;
  sstep = 5;
elseif s > 300,
  step = 10;
  sstep = 2;
elseif s > 100,
  step = 5;
  sstep = 1;
else
  step = 1;
  sstep = 1;
end

kk = 1;

% ----------------------------------------- Put nucleotide numbers outside

for k = 1:s,
  kkk = str2num(File.NT(k).Number);
  if ~isempty(kkk),
    kk = kkk;
  end
  theta = k*tp/s;
  if cos(theta) > 0,
    angle = 180*theta/pi;
    ha    = 'left';
  else
    angle = 180*(theta - pi)/pi;
    ha    = 'right';
  end

  if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1),
    plot([cos(theta) 1.02*cos(theta)], [sin(theta), 1.02*sin(theta)],'k');
    text(1.04*cos(theta), 1.04*sin(theta), File.NT(k).Number,'FontSize',4, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end
  if mod(kk,step) == 0,
    plot([cos(theta) 1.04*cos(theta)], [sin(theta), 1.04*sin(theta)],'k');
    text(1.11*cos(theta), 1.11*sin(theta), File.NT(k).Number, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end

end

% ---------------------------------------- Draw the interactions in color

[y,k] = sort(-abs(c));               % sort by decreasing color, for overlap
i = i(k);
j = j(k);
c = c(k);

for k = 1:length(i),
  if c(k) == 6,
    zUnitCircleArc([cos((i(k)+0.2)*tp/s) cos((j(k)-0.2)*tp/s)], [sin((i(k)+0.2)*tp/s) sin((j(k)-0.2)*tp/s)],c(k),Thickness);
  else
    zUnitCircleArc([cos(i(k)*tp/s) cos(j(k)*tp/s)], [sin(i(k)*tp/s) sin(j(k)*tp/s)],c(k),Thickness);
  end
end
axis equal
axis([-1.2 1.2 -2 1.2]);
axis off

text(-1.2,1.2,File.Filename,'HorizontalAlignment','Left');
text(-1.2,-1.4,'Dark blue chords indicate nested Watson-Crick basepairs');
text(-1.2,-1.5,'Red chords indicate non-nested Watson-Cricky basepairs');
text(-1.2,-1.6,'Cyan chords indicate nested non-Watson-Crick basepairs');
text(-1.2,-1.7,'Green chords indicate non-nested non-Watson-Crick basepairs');
text(-1.2,-1.8,'Yellow chords indicate stacking interactions');
text(-1.2,-1.9,'Magenta chords indicate base-phosphate interactions');

end
