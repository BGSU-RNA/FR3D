% zEnterViewMode(Param,ViewParam) is part of zPairViewer.  It prompts the
% user to specify the view mode, which is coded this way:
% 1  scatter plots
% 2  list of pair parameters
% 3  sorted list of pair parameters
% 4  view individual nucleotide pairs
% 5  view context for these pairs
% When the scatterplot is chosen, it also prompts for the color scheme to
% use to color the points.

function [ViewParam] = zEnterViewMode(Param,ViewParam)

% --------- Set default values for selection parameters ---------------

if ViewParam.FigNum == 2,
  ViewParam.Mode      = 1; 
  ViewParam.Color     = 9;
  ViewParam.Normal    = 0;
  ViewParam.ColorAxis = [-12 30];
  ViewParam.SortKeys  = [];
  ViewParam.Nearby    = 0;
  ViewParam.Sugar     = 0;
  ViewParam.Hydrogen  = 1;
  ViewParam.Sort      = 0;
  ViewParam.az        = 51;
  ViewParam.el        = 14;
  ViewParam.LineStyle = '-';              % default - thick solid lines
  ViewParam.Exemplars = 0;
  ViewParam.ClassLimits = 1;
end

fprintf('\n');

if Param.Expert == 0,
  for k=1:5,
    fprintf('Mode %1d - %s\n',k,zDisplayModeName(k));
  end  
end

if Param.Category == 51,
  ViewParam.Color = 13;
end


fprintf('Enter view mode(s), 9 to sort, 8 for a new pair, 0 to quit [%1d] ', ViewParam.Mode);
inp = input('');
if ~isempty(inp),
  if inp == 9,
    ViewParam = zEnterSortKeys(ViewParam);
    ViewParam.Sort = 1;
  elseif inp == 8,
    ViewParam.Mode = 0;
  elseif inp == 0,
    ViewParam.Mode = 0;
    ViewParam.FigNum = 0;
  else
    ViewParam.Mode = sort(inp);           % do these in a specific order
    ViewParam.Sort = 0;
  end
end

if ((ViewParam.Mode(1) == 1) | (ViewParam.Mode(1) == 6)) & (ViewParam.Sort == 0),          % scatterplot
  if Param.Expert == 0,
    for k=1:13,
      fprintf('Scheme %1d - %s\n',k,zColorSchemeName(k));
    end  
  end

  fprintf('Enter color scheme [%1d] ', ViewParam.Color);
  inp = input('');
  if ~isempty(inp),
    ViewParam.Color = inp;
  end

  fprintf('Enter 1 to show normal vector of second base [%1d] ', ...
           ViewParam.Normal);
  inp = input('');
  if ~isempty(inp),
    ViewParam.Normal = inp;
  end

  fprintf('Enter 1 to show class limits [%1d] ', ...
           ViewParam.ClassLimits);
  inp = input('');
  if ~isempty(inp),
    ViewParam.ClassLimits = inp;
  end

end      
