% zCategoryName(a) returns the category numbered a

function [n] = zCategoryName(a)

switch fix(a),
  case   0, n = 'no classification'
  case   1, n = 'cWW cis Watson-Crick / Watson-Crick';
  case   2, n = 'tWW trans Watson-Crick / Watson-Crick';
  case   3, n = 'cWH cis Watson-Crick / Hoogsteen';
  case  -3, n = 'cHW cis Hoogsteen / Watson-Crick';
  case   4, n = 'tWH trans Watson-Crick / Hoogsteen';
  case  -4, n = 'tHW trans Hoogsteen / Watson-Crick';
  case   5, n = 'cWS cis Watson-Crick / Sugar edge';
  case  -5, n = 'cSW cis Sugar edge / Watson-Crick';
  case   6, n = 'tWS trans Watson-Crick / Sugar edge';
  case  -6, n = 'tSW trans Sugar edge / Watson-Crick';
  case   7, n = 'cHH cis Hoogsteen / Hoogsteen';
  case   8, n = 'tHH trans Hoogsteen / Hoogsteen';
  case   9, n = 'cHS cis Hoogsteen / Sugar edge';
  case  -9, n = 'cSH cis Sugar edge / Hoogsteen';
  case  10, n = 'tHS trans Hoogsteen / Sugar edge';
  case -10, n = 'tSH trans Sugar edge / Hoogsteen';
  case  11, n = 'cSS cis Sugar edge / Sugar edge (second base dominant)';
  case -11, n = 'cSS cis Sugar edge / Sugar edge (first base dominant)';
  case  12, n = 'tSS trans Sugar edge / Sugar edge (second base dominant)';
  case -12, n = 'tSS trans Sugar edge / Sugar edge (first base dominant)';
  case  13, n = 'bif bifurcated cis Watson-Crick / Watson-Crick';
  case -13, n = 'bif bifurcated cis Watson-Crick / Watson-Crick';
  case  14, n = 'mis miscellaneous hand classification';
  case  21, n = 's35 second base faces up, above first';
  case -21, n = 's53 second base faces up, below first';
  case  22, n = 's33 second base faces down, above first';
  case  23, n = 's55 second base faces down, below first';
  case  25, n = 'part of a motif (not implemented)';
  case  26, n = 'sugar stacked on base (not implemented)';
  case  30, n = 'no interaction';
  case  40, n = 'potential stacking above, second base facing up';
  case  41, n = 'potential stacking above, second base facing down';
  case  42, n = 'potential stacking below, second base facing up';
  case  43, n = 'potential stacking below, second base facing down';
  case  44, n = 'potential pairing';
  case  45, n = 'potential pairing, no computer classification';
  case  50, n = 'in given box';
  case  51, n = 'near a pair you specify';
  otherwise, n = 'unknown';
end
