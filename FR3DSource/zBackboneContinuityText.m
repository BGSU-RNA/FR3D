
function [T] = zBackboneContinuityText(c)

switch c,
  case  1, T = 'c35';
  case  2, T = 'c25';
  case -1, T = 'c53';
  case -2, T = 'c52';
  otherwise T = '   ';
end
