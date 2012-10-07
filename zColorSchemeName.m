% zColorSchemeName(ColorScheme) prompts the user to choose a scheme for
% coloring the points plotted by zScatterPairs

function [n] = zColorSchemeName(ColorScheme)

switch ColorScheme,
  case 1, n = 'Blue - cutoff matches; Green - hand matches; Orange - both agree'; 
  case 2, n = 'Cutoff classification -12 to 12';
  case 3, n = 'Cutoff classification 15 to 18';
  case 4, n = 'Cutoff classification -12 to 30';
  case 5, n = 'Pair code'; 
  case 6, n = 'Appropriate angle of rotation';
  case 7, n = 'Gap'; 
  case 8, n = 'Blue - cutoff matches; Green - exemplar matches; Orange - both agree';
  case 9, n = 'Distance to nearest exemplar';
  case 10, n = 'Exemplar classification';
  case 11, n = 'Degree of stacking overlap';
  case 12, n = 'Subcategory of first category';
  case 13, n = 'Distance to specified pair';
end
