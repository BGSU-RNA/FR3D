% xScatterPairs displays multiple 3d scatter plots of pair parameters, with a menu that allows one to color points in different orders, or display different nucleotides from the candidates

function [ppp] = xScatterPairs(Search,N1,N2,ViewParam,Param)

warning off

if nargin < 3,
  N1 = 1;
  N2 = 2;
end

% --------------------------------------------- 
% --------------------------------------------- Collect basic data

Candidates = Search.Candidates;
[L,t] = size(Candidates);
N     = t - 1;                                     % number of nucleotides
File  = Search.File;
CL    = zClassLimits;                              % read ClassLimits matrix

ppp = 1:L;                                         % default ordering

if exist('PairExemplars.mat','file') > 0,
  load('PairExemplars','Exemplar');
else
  Exemplar = [];
end

% --------------------------------------------- Provide memory space

for i = 1:N,
  for j = 1:N,
    SavedPairs(i,j).Classified = 0;
  end
end

% --------------------------------------------- Set defaults

if nargin < 4,
  ViewParam.Color  = 6;
  ViewParam.FigNum = 1;
  ViewParam.Normal = 0;
  ViewParam.ClassLimits = 1;
end

if nargin < 5,
  Param.Decimal    = 1;
end

%---------------------------------------------- Loop for the menu

Stop       = 0;
Recolor    = 1;
Reclassify = 1;
Replot     = 1;

while Stop == 0,


% --------------------------------------------- Analyze pairs of nucleotides

if Reclassify > 0,

 if SavedPairs(N1,N2).Classified == 0,

  pc = 1;                                          % pair counter
  FirstBaseCodes = [];
  Paircode       = [];
  fprintf('Nucleotide %3d at the origin, %3d away\n', N1, N2);

  for k = 1:L,                                     % Loop through candidates
    f  = Candidates(k,N+1);                        % file number
    i1 = Candidates(k,N1);                         % first nucleotide index
    i2 = Candidates(k,N2);                         % second nucleotide index  

    Ni = File(f).NT(i1);
    Nj = File(f).NT(i2);

    [p,s] = zClassifyPair(Ni,Nj,CL,Exemplar,1);

    p.NT1 = Ni;
    p.NT2 = Nj;
    p.Filename = File(f).Filename;
    p.C1pC1p   = zDistance(Ni.Sugar(1,:),Nj.Sugar(1,:));  % C1' - C1' distance
    p.Resol    = File(f).Info.Resolution;
    if isempty(p.Resol)
      p.Resol = NaN;
    end

    if ~isempty(p),
      Pair(pc) = p;                                % store this pair
      pc = pc + 1;                                 % increment pair counter
      FirstBaseCodes(Ni.Code) = 1;                 % this base was seen
      Paircode(p.Paircode)    = 1;                 % this paircode was seen
    end

  end

  modcode = find(FirstBaseCodes);                  % codes of first bases seen
  Paircodes = find(Paircode);

  SavedPairs(N1,N2).Pair = Pair;
  SavedPairs(N1,N2).modcode = modcode;
  SavedPairs(N1,N2).Paircodes = Paircodes;

 else

  Pair      = SavedPairs(N1,N2).Pair;
  modcode   = SavedPairs(N1,N2).modcode;
  Paircodes = SavedPairs(N1,N2).Paircodes;
  
 end

end

%---------------------------------------------- Set the colors for the points

if Recolor > 0,

 for k = 1:length(Pair),                                   % Loop through pairs
  p = Pair(k);

  % -------- Round categories to nearest integer, if desired ------------- 

  if Param.Decimal == 1,                  % select with decimal places in categ
    pClass = p.Edge;
    eClass = p.Classes(1);
  else
    pClass = fix(p.Edge);
    eClass = fix(p.Classes(1));
  end

  switch ViewParam.Color,
    case  1, Color(k) = k;
    case {2, 3, 4}, Color(k) = p.Edge;
    case  5, Color(k) = p.Paircode;
    case  6, Color(k) = p.Ang;
    case  7, Color(k) = p.Gap;
    case  8, Color(k) = p.Normal(3);
    case  9, Color(k) = p.Distances(1);
    case 11, Color(k) = p.StackingOverlap;
    case 12, Color(k) = p.Edge;
    case 13, Color(k) = p.Displ(1);
    case 14, Color(k) = p.Displ(2);
    case 15, Color(k) = p.Displ(3);
    otherwise, Color(k) = 1;
  end

 end

  [y,i] = sort(Color);
  ppp   = i;                                    % permutation
  Pair  = Pair(i);
  Color = Color(i);

end

% --------------------------------------------- Set color axis

switch ViewParam.Color,
  case 1, ColorAxis =  [min(Color) max(Color)];
  case 2, ColorAxis =  [1 12];
  case 3, ColorAxis =  [15 19];
  case 4, ColorAxis =  [-12 30];
  case 5, ColorAxis =  [1 16];
  case 6, ColorAxis =  [min(Color) max(Color)];      % color by angle
  case 7, ColorAxis =  [min(Color) max(Color)];
  case 8, ColorAxis =  [-1 3];
  case 9, ColorAxis =  [0 4];
  case 10, ColorAxis = [-12 12];
  case 11, ColorAxis = [0 10];
  case 12, ColorAxis = [min(Color) max(Color)];
  case 13, ColorAxis =  [min(Color) max(Color)];
  case 14, ColorAxis =  [min(Color) max(Color)];
  case 15, ColorAxis =  [min(Color) max(Color)];
  otherwise, ColorAxis =  [0 10];
end

%---------------------------------------------- Plot displacements

if Replot > 0,

figure(ViewParam.FigNum)
clf
figure(ViewParam.FigNum)
clf

%set(gcf,'Renderer','OpenGL');     % fast rotation
%set(gcf,'Renderer','zbuffer')
%set(gcf,'Renderer','painters')

for k = 1:length(Pair),                              % Loop through pairs
  p      = Pair(k);
  c(k)   = Color(k);
  e(k,:) = p.Displ;

  if ViewParam.Normal == 1,
    v = p.Normal/3;                                      % add normal vector
    plot3([e(k,1) e(k,1)+v(1)], [e(k,2) e(k,2)+v(2)], [e(k,3) e(k,3)+v(3)], 'b');
    hold on
  end
end

scatter3(e(:,1),e(:,2),e(:,3),18*ones(size(c)),c,'filled')

if max(modcode) == min(modcode),             % all have same first base
  zPlotStandardBase(modcode);                % plot base at the origin
  Lett = 'ACGU';  
  Title =[Lett(modcode) ' shown at the origin, N1/N9 atom of second base shown by dots'];
else
  Title = ['N1/N9 atom of second base shown by dots'];
end

% ------------------ Display class limits, if desired

if isfield(ViewParam,'ClassLimits'),
 if ViewParam.ClassLimits == 1,                  % show all boxes
   B = CL(:,:,Paircodes(1));                     % use limits for this paircode

   for row = 1:length(B(:,1)),
     hold on
     if (abs(B(row,1)) < 14) & (abs(B(row,1)) > 0),
       zSquare([B(row,[2 4]) 0],[B(row,[3 5]) 0],'k');
       tex = zEdgeText(B(row,1),1,Paircodes(1));
       tex = tex(1:4);
       text(B(row,2)+0.1,B(row,4)+0.2,0,tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
     end
   end
 elseif ViewParam.ClassLimits == 2,              % just show one box
   ClassLimits = zClassLimits;                   % load class limits
   B = CL(:,:,Paircodes(1));       % use limits for this paircode

   for row = 1:length(B(:,1)),
     hold on
     if any(fix(B(row,1)) == fix(Param.Category)),
       zSquare([B(row,[2 4]) 0],[B(row,[3 5]) 0],'k');
       tex = zEdgeText(B(row,1),1,Paircodes(1));
       tex = tex(1:4);
       text(B(row,2)+0.1,B(row,4)+0.2,0,tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
     end
   end
 end
end

zoom on

caxis(ColorAxis);
view(2)
axis equal

if ViewParam.Normal == 1,
 Title = [Title ', lines are the normal vector of second base'];
end
title(Title);
xlabel('Perpendicular to glycosidic bond');
ylabel('Parallel to glycosidic bond');
zlabel('Vertical with respect to first base');

%----------------------------------------------------- Plot angle and normal
figure(ViewParam.FigNum + 1)
clf
figure(ViewParam.FigNum + 1)
clf

%set(gcf,'Renderer','OpenGL');
%set(gcf,'Renderer','zbuffer')

cut = 90;                                             % cut at angle -cut

for k = 1:length(Pair),                               % Loop through pairs
  p = Pair(k);
  e = p.Displ;

  c = Color(k);

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
    B = CL(:,:,Paircodes(1));     % use limits for this paircode

    B(:,10:11) = mod(B(:,10:11)+cut,360)-cut;
    B(:,8:9)   = B(:,8:9) + 0.1*(rand(size(B(:,8:9)))-0.5).*(abs(B(:,8:9))>1);

    for row = 1:length(B(:,1)),
      hold on
      if (abs(B(row,1)) < 14) & (abs(B(row,1)) > 0),
        tex = zEdgeText(B(row,1),1,Paircodes(1));
        tex = tex(1:4);
        if (B(row,10) < B(row,11)),
          zSquare([B(row,[10 8]) 0],[B(row,[11 9]) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
        else
          zSquare([B(row,[10 8]) 0],[270 B(row,9) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
          zSquare([-90 B(row,8) 0],[B(row,[11 9]) 0],'k');
          text(-86,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
        end
      end
    end
  elseif ViewParam.ClassLimits == 2,            % show only one box
    B = CL(:,:,Paircodes(1));     % use limits for this paircode
    B(:,10:11) = mod(B(:,10:11)+cut,360)-cut;
    B(:,8:9)   = B(:,8:9) + 0.1*(rand(size(B(:,8:9)))-0.5).*(abs(B(:,8:9))>1);

    for row = 1:length(B(:,1)),
      hold on
      if any(fix(B(row,1)) == fix(Param.Category)),
        tex = zEdgeText(B(row,1),1,Paircodes(1));
        tex = tex(1:4);
        if (B(row,10) < B(row,11)),
          zSquare([B(row,[10 8]) 0],[B(row,[11 9]) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
        else
          zSquare([B(row,[10 8]) 0],[270 B(row,9) 0],'k');
          text(B(row,10)+4,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
          zSquare([-90 B(row,8) 0],[B(row,[11 9]) 0],'k');
          text(-86,B(row,8)+0.08,B(row,11),tex,'horizontalalignment','left','fontweight','bold','FontSize',8);
        end
      end
    end
    v = axis;
    axis([-cut 360-cut -1.1 1.1 -2 2]);
  end
end

zoom on

%----------------------------------------------------- Plot angle and normal

if ViewParam.Color > 1,
  figure(ViewParam.FigNum + 2)
  clf
  figure(ViewParam.FigNum + 2)
  clf

  hist(Color,30)
  switch ViewParam.Color,
    case {2, 3, 4}, title('Histogram of interaction code');
    case  5, title('Histogram of paircode'); % Color(k) = p.Paircode;
    case  6, title('Histogram of angle of rotation'); % Color(k) = p.Ang;
    case  7, title('Histogram of gap'); % Color(k) = p.Gap;
    case  8, title('Histogram of third component of normal vector');
    case  9, title('Histogram of distance'); % Color(k) = p.Distances(1);
    case 11, title('Histogram of stacking overlap'); % Color(k) = p.StackingOverlap;
    case 12, title('Histogram of interaction code'); % Color(k) = p.Edge;
    case 13, title('Histogram of perpendicual displacement'); % Color(k) = p.Displ(1);
    case 14, title('Histogram of parallel displacement'); % Color(k) = p.Displ(2);
    case 15, title('Histogram of vertical displacement'); % Color(k) = p.Displ(3);
  end
end


end
% ---------------------------------------------------- Display the menu

  k=menu('Scatterplot controls',...
         'Increment first nucleotide',...             % 1
         'Decrement first nucleotide',...             % 2
         'Increment second nucleotide',...            % 3
         'Decrement second nucleotide',...            % 4
         'Color by angle',...                         % 5
         'Color by pair code',...                     % 6
         'Color by interaction',...                   % 7
         'Color by subcategory',...                   % 8
         'Color by displacement perp to bond',...     % 9
         'Color by displacement parallel to bond',... % 10
         'Color by vertical displacement',...         % 11
         'Color by position in list',...              % 12
         'Color by gap',...                           % 13
         'Color by third component of normal',...     % 14
         'Toggle normal vector',...                   % 15
         'Toggle category limits',...                 % 16
         'List pair parameters',...                   % 17
         'Quit');                                     % 18

  switch k,
  case 1, N1 = N1 + 1;
          if N1 > N,
            N1 = 1;
          end
          if N1 == N2,
            N1 = N2 + 1;
          end
          if N1 > N,
            N1 = 1;
          end
  case 2, N1 = N1 - 1;
          if N1 < 1,
            N1 = N;
          end
          if N1 == N2,
            N1 = N2 - 1;
          end
          if N1 < 1,
            N1 = N;
          end
  case 3, N2 = N2 + 1;
          if N2 > N,
            N2 = 1;
          end
          if N2 == N1,
            N2 = N2 + 1;
          end
          if N2 > N,
            N2 = 1;
          end
  case 4, N2 = N2 - 1;
          if N2 < 1,
            N2 = N;
          end
          if N2 == N1,
            N2 = N2 - 1;
          end
          if N2 < 1,
            N2 = N;
          end
  case  5, ViewParam.Color = 6;
  case  6, ViewParam.Color = 5;
  case  7, ViewParam.Color = 2;
  case  8, ViewParam.Color = 12;
  case  9, ViewParam.Color = 13;
  case 10, ViewParam.Color = 14;
  case 11, ViewParam.Color = 15;
  case 12, ViewParam.Color = 1;
  case 13, ViewParam.Color = 7;
  case 14, ViewParam.Color = 8;
  case 15, ViewParam.Normal = 1 - ViewParam.Normal;
           Recolor = 0;
  case 16, ViewParam.ClassLimits = 1 - ViewParam.ClassLimits;
           Recolor = 0;
  case 17, xListPairs(Pair,ViewParam);
           Replot = 0;
  case 18, Stop = 1;
  end

  if (k >= 1) && (k <= 4) && (N == 2),     % switch order
    N3 = N1;
    N1 = N2;
    N2 = N3;
  end

  if k < 17,
    Replot  = 1;
  end

  if k < 13,
    Recolor  = 1;
  end

  if k < 5,
    Reclassify = 1;
  end

end

FigsDone = 2;

warning on