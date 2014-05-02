% zDisplayAA(File,AAList,ViewParam) plots amino acids.  It can be called in several ways, for example,
% zDisplayAA('1s72',{'27','28','29'},ViewParam), and it will load the
% named datafile and plot the amino acids by amino acid number.  Or, if
% data files have already been loaded, one can use zDisplayAA(File(1),[32
% 34 35],ViewParam) to plot the amino acids in File(1) having indices 32,
% 34, and 35.  Defaults for ViewParam are defined in zDisplayAA; see there
% for the fields of ViewParam.
% One can also use ranges of amino acid numbers, as in
% zDisplayAA('rr0033_23S',{'2548:2555','2557','2559:2566'},VP);

function [void] = zDisplayAA(File,AAList,ViewParam)

% set default values

VP.Sugar     = 0;
VP.LabelSugar= 0;
VP.az        = 51;
VP.el        = 14;
VP.LineStyle = '-';              % default - thick solid lines
VP.LineThickness = 2;          
VP.AtOrigin  = 0;                % rotate all so first is at origin
VP.AATitle     = 0;                % title with amino acid numbers, filename
VP.Grid      = 1;                % add a grid to the graph
VP.FontSize  = 10;               % will use Matlab's default unless overridden
VP.Rotation  = eye(3);
VP.Shift     = zeros(1,3);
VP.LabelBases = 10;              % font size; use 0 to not label bases at all
VP.GlycoAtomSize = 28;           % size to mark the glycosidic atom
VP.WritePDB      = 0;
VP.Animate   = 0;
VP.AnimationFilename = '';
VP.ShowBeta  = 0;                % show beta factor for each atom
VP.AABackboneTrace = 0;

if nargin == 1,
  ViewParam = VP;
end

if nargin == 2,
  ViewParam = VP;
end

% replace defaults with defined values

if isfield(ViewParam,'Sugar'),
  VP.Sugar = ViewParam.Sugar;
end

if isfield(ViewParam,'LabelSugar'),
  VP.LabelSugar = ViewParam.LabelSugar;
end

if isfield(ViewParam,'az'),
  VP.az = ViewParam.az;
end

if isfield(ViewParam,'el'),
  VP.el = ViewParam.el;
end

if isfield(ViewParam,'LineStyle'),
  VP.LineStyle = ViewParam.LineStyle;
end

if isfield(ViewParam,'LineThickness'),
  VP.LineThickness = ViewParam.LineThickness;
end

if isfield(ViewParam,'Color'),
  VP.Color = ViewParam.Color;
end

if isfield(ViewParam,'AtOrigin'),
  VP.AtOrigin = ViewParam.AtOrigin;
end

if isfield(ViewParam,'AATitle'),
  VP.AATitle = ViewParam.AATitle;
end

if isfield(ViewParam,'Grid'),
  VP.Grid = ViewParam.Grid;
end

if isfield(ViewParam,'FontSize'),
  h=gca;
  set(h,'FontSize',ViewParam.FontSize);
end

if isfield(ViewParam,'Rotation'),
  VP.Rotation = ViewParam.Rotation;
end

if isfield(ViewParam,'Shift'),
  VP.Shift = ViewParam.Shift;
end

if isfield(ViewParam,'LabelBases'),
  VP.LabelBases = ViewParam.LabelBases;
end

if isfield(ViewParam,'GlycoAtomSize'),
  VP.GlycoAtomSize = ViewParam.GlycoAtomSize;
end

if isfield(ViewParam,'WritePDB'),
  VP.WritePDB = ViewParam.WritePDB;
end

if isfield(ViewParam,'Animate'),
  VP.Animate = ViewParam.Animate;
end

if isfield(ViewParam,'AnimationFilename'),
  VP.AnimationFilename = ViewParam.AnimationFilename;
end

if isfield(ViewParam,'Colors'),
  VP.Colors = ViewParam.Colors;
end

if isfield(ViewParam,'ShowBeta'),
  VP.ShowBeta = ViewParam.ShowBeta;
end

if isfield(ViewParam,'AABackboneTrace'),
  VP.AABackboneTrace = ViewParam.AABackboneTrace;
end

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1 || isempty(AAList),
  AAList = 1:length(File.AA);                  % display them all
end

% if AAList is a cell array of numbers, look up the indices

if strcmp(class(AAList),'char'),
  AAList = {AAList};
end

if strcmp(class(AAList),'cell'),
  Indices = zIndexLookup(File,AAList);
else
  Indices = AAList;
end

% plot the amino acids, at the origin or in original positions

%set(gcf,'Renderer','OpenGL');

if VP.AtOrigin > 0,
  v = fix(VP.AtOrigin);
  R = File.AA(Indices(v)).Rot;             % Rotation matrix for first base
  S = File.AA(Indices(v)).Loc(1,:);        % Location of glycosidic atom
else
  R = VP.Rotation;
  S = VP.Shift;
end

for j=1:length(Indices),                 % Loop through all amino acids
  k = Indices(j);                        % index of current amino acid
  if ~isempty(File.AA(k).Loc),
    zPlotOneAARotated(File.AA(k),VP,R,S);  % connect amino acids

    tcol = [ 205   41    144]/205;             % default
    tcol = [1 1 0.5];

    if VP.AABackboneTrace > 0 && any(k+1 == Indices) && ~isempty(File.AA(k+1).Loc),
      p1 = (File.AA(k).Loc(1,:) - S) * R;
      p2 = (File.AA(k+1).Loc(1,:) - S) * R;
      if norm(p1-p2) < 8,
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)],'color',tcol,'linewidth',5);
      end
    end

    col = [ 205   41    144]/255;             % default

    if any(k+1 == Indices) && ~isempty(File.AA(k+1).Loc),  % connect protein backbone
      p1 = (File.AA(k).Loc(4,:) - S) * R;
      p2 = (File.AA(k+1).Loc(1,:) - S) * R;
      if norm(p1-p2) < 3,
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)],'color',col,'linewidth',2);
      end
    end
  end
end

Title = strcat(File.AA(Indices(1)).Unit,File.AA(Indices(1)).Number);
for j=2:length(Indices),
  nt = File.AA(Indices(j));
  Title = strcat(Title,'-',nt.Unit,nt.Number);
end;

if isfield(File,'Filename'),
  FN = File.Filename;
else
  FN = '';
end

Title = strcat(Title,[' ' strrep(FN,'_','\_')]);

if VP.AATitle > 0,
  title(Title);
end
axis equal
view(VP.az,VP.el);

if VP.Grid == 1,
  grid on
else
  grid off
end

rotate3d on

if VP.WritePDB == 1,
  if VP.AtOrigin > 0,
    v = fix(VP.AtOrigin);
    G = eye(3);
    R = File.AA(Indices(v)).Rot * G(:,[2 3 1]);% Rotation matrix for first base
    S = File.AA(Indices(v)).Loc(1,:);        % Location of glycosidic atom
  else
    R = VP.Rotation;
    S = VP.Shift;
  end

  a = 1;
  fid = fopen('PDBFile.pdb','w');

  for j=1:length(Indices),                 % Loop through all amino acids
    k = Indices(j);                        % index of current amino acid
%    a = zWriteAAPDB(fid,File.AA(k),a,0,R,S);  % not written yet!
  end

  fclose(fid);
end

if VP.Animate == 1,
  NumFrames = 90;
  S = zeros(1,3);
  for j = 1:length(Indices);
    S = S + File.AA(Indices(j)).Center;
  end
  S = S / length(Indices);

  th = 0.1;
  Spin = [[cos(th) sin(th) 0]; [-sin(th) cos(th) 0]; [0 0 1]];

  if isfield(ViewParam,'Rotation'),
    R = VP.Rotation;
  else
    R = File.AA(Indices(1)).Rot * Spin([3 1 2], [3 1 2]);
  end

  maxradius = 0;
  minheight = 0;
  maxheight = 0;
  for j = 1:length(Indices);
    L  = length(File.AA(Indices(j)).Loc(:,1));         % Number of base atoms
    V  = (File.AA(Indices(j)).Loc - ones(L,1)*S) * R;
    maxradius = max(maxradius,max(sqrt(V(:,1).^2+V(:,2).^2)));
    minheight = min(minheight,min(V(:,3)));
    maxheight = max(maxheight,max(V(:,3)));
  end
  maxradius = maxradius + 1;
  minheight = minheight - 1;
  maxheight = maxheight + 1;

  th = 2*pi/NumFrames;
  Spin = [[cos(th) sin(th) 0]; [-sin(th) cos(th) 0]; [0 0 1]];

  for g = 1:NumFrames,

    R = R * Spin;
    clf
    for j=1:length(Indices),                 % Loop through all amino acids
      k = Indices(j);                        % index of current amino acid
      if isfield(VP,'Colors'),
        VP.Color = VP.Colors(k,:);
      end
      zPlotOneAARotated(File.AA(k),VP,R,S);
    end
  
    Title = strcat(File.AA(Indices(1)).Unit,File.AA(Indices(1)).Number);
    for j=2:length(Indices),
      nt = File.AA(Indices(j));
      Title = strcat(Title,'-',nt.Unit,nt.Number);
    end;

    if isfield(File,'Filename'),
      FN = File.Filename;
    else
      FN = '';
    end
  
    Title = strcat(Title,[' ' strrep(FN,'_','\_')]);
  
    axis equal
    axis([-maxradius maxradius -maxradius maxradius minheight maxheight]);
    set(gca,'TickLength',[0 0]);
    grid off  
%    axis off
    view(0,0);
    drawnow

    Image = getframe;		% get current pseudo color image
    Image.cdata = Image.cdata(1:(end-1),2:end,:);
    P = frame2im(Image);		% Convert to a image representation that Matlab can handle
    filename = ['Animation/' VP.AnimationFilename num2str(g+100) '.bmp'];
    imwrite(P,filename, 'bmp')  % Finally write the individual images

  end
end
