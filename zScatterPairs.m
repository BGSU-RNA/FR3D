% zScatterPairs displays multiple 3d scatter plots of pair parameters

function [FigsDone] = zScatterPairs(File,SP,Param,ViewParam)

T = cat(1,SP(:).Color);

switch ViewParam.Color,
  case 1, ColorAxis =  [0 8];
          fprintf('Blue - computer matches; Green - hand matches; Orange - both match\n');
  case 2, ColorAxis =  [1 12];
%  case 2, ColorAxis =  [-12 12];
  case 3, ColorAxis =  [15 19];
  case 4, ColorAxis =  [-12 30];
  case 5, ColorAxis =  [1 16];
  case 8, ColorAxis =  [0 8];
          fprintf('Blue - cutoff matches; Green - exemplar matches; Orange - both match\n');
  case 9, ColorAxis =  [0 4];
  case 10, ColorAxis = [-12 12];
  case 11, ColorAxis = [0 10];
  case 12, ColorAxis = [Param.Category(1) Param.Category(1)+0.3];
  otherwise, ColorAxis =  [min(T) max(T)+0.01];
end

%-----------------------------------------------------------------------
figure(ViewParam.FigNum)
clf

%set(gcf,'Renderer','OpenGL');
set(gcf,'Renderer','zbuffer')

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  e = p.Displ;

  scatter3(e(1),e(2),e(3),18,c,'filled')
  hold on

  if ViewParam.Normal == 1,
    v = p.Normal/3;                                      % add normal vector
    plot3([e(1) e(1)+v(1)], [e(2) e(2)+v(2)], [e(3) e(3)+v(3)], 'b');
  end
end

modcode = mod(Param.Paircode,4) + 4*(mod(Param.Paircode,4)==0);

if max(modcode) == min(modcode),             % all have same first base
  zPlotStandardBase(modcode);
  Lett = 'ACGU';  
  Title =[Lett(modcode) ' shown at the origin, N1/N9 atom of second base shown by dots'];
else
  Title = ['N1/N9 atom of second base shown by dots'];
end

if isfield(ViewParam,'ClassLimits'),
 if ViewParam.ClassLimits == 1,                  % show all boxes
   ClassLimits = zClassLimits;                   % load class limits
   B = ClassLimits(:,:,Param.Paircode(1));       % use limits for this paircode

   for row = 1:length(B(:,1)),
     hold on
     if (abs(B(row,1)) < 14) & (abs(B(row,1)) > 0),
       zSquare([B(row,[2 4]) 0],[B(row,[3 5]) 0],'k');
%       text(B(row,2)+0.2,B(row,4)+0.35,B(row,7)+0.2,num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
       text(B(row,2)+0.2,B(row,4)+0.35,0,num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
     end
   end
 elseif ViewParam.ClassLimits == 2,              % just show one box
   ClassLimits = zClassLimits;                   % load class limits
   B = ClassLimits(:,:,Param.Paircode(1));       % use limits for this paircode

   for row = 1:length(B(:,1)),
     hold on
     if any(fix(B(row,1)) == fix(Param.Category)),
       zSquare([B(row,[2 4]) 0],[B(row,[3 5]) 0],'k');
%       text(B(row,2)+0.2,B(row,4)+0.35,B(row,7)+0.2,num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
       text(B(row,2)+0.2,B(row,4)+0.35,0,num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
     end
   end
 end
end

zoom on

caxis(ColorAxis);
view(2)
axis square

if ViewParam.Normal == 1,
 Title = [Title ', lines are the normal vector of second base'];
end
title(Title);
xlabel('Perpendicular to glycosidic bond');
ylabel('Parallel to glycosidic bond');
zlabel('Vertical with respect to first base');

%----------------------------------------------------------------------
figure(ViewParam.FigNum + 1)
clf

cut = 90;

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  scatter3(mod(p.Ang+cut,360)-cut,p.Normal(3),p.Gap,18,c,'filled')
  hold on
end

xlabel('Angle of rotation (in degrees)');
ylabel('Vertical component of normal');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

view(2)

if isfield(ViewParam,'ClassLimits'),
  if ViewParam.ClassLimits == 1,
    ClassLimits = zClassLimits;                 % load class limits
    B = ClassLimits(:,:,Param.Paircode(1));     % use limits for this paircode

    B(:,10:11) = mod(B(:,10:11)+cut,360)-cut;
    B(:,8:9)   = B(:,8:9) + 0.1*(rand(size(B(:,8:9)))-0.5).*(abs(B(:,8:9))>1);

    for row = 1:length(B(:,1)),
      hold on
      if (abs(B(row,1)) < 14) & (abs(B(row,1)) > 0),
        if (B(row,10) < B(row,11)),
          zSquare([B(row,[10 8]) 0],[B(row,[11 9]) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
        else
          zSquare([B(row,[10 8]) 0],[270 B(row,9) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
          zSquare([-90 B(row,8) 0],[B(row,[11 9]) 0],'k');
          text(-86,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
        end
      end
    end
  elseif ViewParam.ClassLimits == 2,            % show only one box
    ClassLimits = zClassLimits;                 % load class limits
    B = ClassLimits(:,:,Param.Paircode(1));     % use limits for this paircode
    B(:,10:11) = mod(B(:,10:11)+cut,360)-cut;
    B(:,8:9)   = B(:,8:9) + 0.1*(rand(size(B(:,8:9)))-0.5).*(abs(B(:,8:9))>1);

    for row = 1:length(B(:,1)),
      hold on
      if any(fix(B(row,1)) == fix(Param.Category)),
        if (B(row,10) < B(row,11)),
          zSquare([B(row,[10 8]) 0],[B(row,[11 9]) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
        else
          zSquare([B(row,[10 8]) 0],[270 B(row,9) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
          zSquare([-90 B(row,8) 0],[B(row,[11 9]) 0],'k');
          text(-86,B(row,8)+0.08,B(row,11),num2str(B(row,1)),'horizontalalignment','left','fontweight','bold','FontSize',15);
        end
      end
    end
    v = axis;
    axis([-cut 360-cut -1.1 1.1 -2 2]);
  end
end

zoom on

FigsDone = 2;
return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 2)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  h = p.Hydrogen;
  if length(h) > 0
    hm = min(cat(1,h(:).Distance));
  else
    hm = 0;
  end
%  scatter3(SP(k).MinDist,p.Displ(3),p.Gap,18,c,'filled')
  scatter3(hm,p.Displ(3),p.Gap,18,c,'filled')
  hold on
end

xlabel('Minimum distance');
ylabel('Vertical displacement');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

FigsDone = 3;

return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 2)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  if length(p.Hydrogen) > 0,
    hmin = min(cat(1,p.Hydrogen(:).Angle));
    hmax = max(cat(1,p.Hydrogen(:).Angle));
    scatter3(hmin,p.Gap,hmax,18,c,'filled')
  end
  hold on
end

xlabel('Minimum hydrogen angle');
ylabel('Gap');
zlabel('Maximum hydrogen angle');
title('');
caxis(ColorAxis);
grid on

FigsDone = 3;

return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 3)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  scatter3(SP(k).C1pC1p,p.Ang,p.Gap,18,c,'filled')
  hold on
end

xlabel('C1*-C1* distance');
ylabel('Angle of rotation');
zlabel('Gap');
title('');
caxis(ColorAxis);
grid on

FigsDone = 4;

return

%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 4)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  L = length(p.Hydrogen); 
  if L > 0,
    a2 = p.Hydrogen(min(2,L)).Angle;
    a3 = p.Hydrogen(min(3,L)).Angle;
    scatter3(p.Hydrogen(1).Angle,a2,a3,18,c,'filled')
    hold on
  end
end

xlabel('First hydrogen bond');
ylabel('Second hydrogen bond, if any, otherwise first');
zlabel('Third hydrogen bond, if any, otherwise first');
title('');
caxis(ColorAxis);
grid on

FigsDone = 5;


%-----------------------------------------------------------------------
figure(ViewParam.FigNum + 5)
clf

for k = 1:length(SP),                               % Loop through pairs
  p = File(SP(k).Filenum).Pair(SP(k).PairIndex);    % Current pair
  c = SP(k).Color;
  L = length(p.Hydrogen); 
  if L > 0,
    a2 = p.Hydrogen(min(2,L)).Distance;
    a3 = p.Hydrogen(min(3,L)).Distance;
    scatter3(p.Hydrogen(1).Distance,a2,a3,18,c,'filled')
    hold on
  end
end

xlabel('First hydrogen bond length');
ylabel('Second hydrogen bond length, if any, otherwise first');
zlabel('Third hydrogen bond length, if any, otherwise first');
title('');
caxis(ColorAxis);
grid on

FigsDone = 6;


return

%---------------------------------------------------------------------

scatter3(RotAx(G,1),RotAx(G,2),RotAng(G),18,T,'filled')
title('Two components of axis of rotation, vertical is angle')
xlabel('Component 1 of axis of rotation')
ylabel('Component 2 of axis of rotation')
zlabel('Angle of rotation')
axis([-1 1 -1 1 -90 270]);
caxis(ColorAx);

figure(4)
clf
scatter3(MinDist(G),PlaneAng(G),Ang(G),18,T,'filled')
xlabel('Minimum distance');
ylabel('Angle between planes');
zlabel('Appropriate angle of rotation');
caxis(ColorAx);

