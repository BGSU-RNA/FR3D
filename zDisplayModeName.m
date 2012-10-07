% zDisplayModeName(DisplayMode) returns a text string telling the name of
% the various display modes in zPairViewer

function [n] = zDisplayModeName(DisplayMode)

switch DisplayMode,
  case 1, n = 'Scatterplots of pairs'; 
  case 2, n = 'List of pair parameters'; 
  case 3, n = 'Pair parameter statistics'; 
  case 4, n = 'View individual pairs and hand classify'; 
  case 5, n = 'View context for these pairs'; 
end
