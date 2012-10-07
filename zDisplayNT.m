% zDisplayNT(File,NTList,ViewParam) is a general-purpose nucleotide plotting program.  It can be called in several ways, for example,
% zDisplayNT('1s72',{'27','28','29'},ViewParam), and it will load the
% named datafile and plot the nucleotides by nucleotide number.  Or, if
% data files have already been loaded, one can use zDisplayNT(File(1),[32
% 34 35],ViewParam) to plot the nucleotides in File(1) having indices 32,
% 34, and 35.  Defaults for ViewParam are defined in zDisplayNT; see there
% for the fields of ViewParam.
% One can also use ranges of nucleotide numbers, as in
% zDisplayNT('rr0033_23S',{'2548:2555','2557','2559:2566'},VP);

function [void] = zDisplayNT(File,NTList,ViewParam)

% set default values

VP.Sugar     = 0;
VP.az        = 51;
VP.el        = 14;
VP.LineStyle = '-';              % default - thick solid lines
VP.LineThickness = '2';          
VP.AtOrigin  = 0;                % rotate all so first is at origin
VP.Title     = 1;                % title with nucleotide numbers, filename
VP.Grid      = 1;                % add a grid to the graph
VP.FontSize  = 10;               % will use Matlab's default unless overridden
VP.Rotation  = eye(3);
VP.Shift     = zeros(1,3);
VP.LabelBases = 10;

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

if isfield(ViewParam,'ConnectSugar'),
  VP.ConnectSugar = ViewParam.ConnectSugar;
else
  if VP.Sugar > 0,
    VP.ConnectSugar = 1;
  else
    VP.ConnectSugar = 0;
  end
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

if isfield(ViewParam,'Title'),
  VP.Title = ViewParam.Title;
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

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1,
  NTList = 1:File.NumNT;                  % display them all
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

% plot the nucleotides, at the origin or in original positions

%set(gcf,'Renderer','OpenGL');

if VP.AtOrigin == 1,
  R = File.NT(Indices(1)).Rot;             % Rotation matrix for first base
  S = File.NT(Indices(1)).Fit(1,:);        % Location of glycosidic atom
else
  R = VP.Rotation;
  S = VP.Shift;
end

for j=1:length(Indices),                 % Loop through all nucleotides
  k = Indices(j);                        % index of current nucleotide
  zPlotOneNTRotated(File.NT(k),VP,R,S);
end

% connect the sugars, if sugars are being plotted

if (VP.Sugar > 0) & (VP.ConnectSugar > 0),
  [D,i] = sort(Indices);
  for j=1:(length(Indices)-1),
    if D(j+1)-D(j) == 1,
      A = [File.NT(D(j)).Sugar(5,:); File.NT(D(j+1)).Sugar(10,:)];
      AA = (A - ones(2,1)*S) * R;
      plot3(AA(:,1), AA(:,2), AA(:,3),'k','LineWidth',2,'LineStyle','-');
    end
  end
end
  
Title = strcat(File.NT(Indices(1)).Base,File.NT(Indices(1)).Number);
for j=2:length(Indices),
  nt = File.NT(Indices(j));
  Title = strcat(Title,'-',nt.Base,nt.Number);
end;

Title = strcat(Title,[' ' strrep(File.Filename,'_','\_')]);

title(Title);
axis equal
view(VP.az,VP.el);

if VP.Grid == 1,
  grid on
else
  grid off
end

rotate3d on
