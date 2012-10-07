
function [T] = zBackboneText(bcc)

if bcc > 0,
  zBackboneCodes
  T = Codes{fix(bcc)};
else
  T = ' -';
end
