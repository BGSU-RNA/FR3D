function [n] = zGroupNames(group)

switch group,
  case 1, n = 'Computer classification matches category';
  case 2, n = 'Hand classification matches category';
  case 3, n = 'Either cutoff or hand classification matches'; 
  case 4, n = 'Computer classification matches but hand differs'; 
  case 5, n = 'Hand matches but computer differs';
  case 6, n = 'Computer classification matches and pair is sequential';
  case 7, n = 'Nearest exemplar matches';
  case 8, n = 'Nearest exemplar matches but cutoff classification differs';
  case 9, n = 'Computer classification matches but nearest exemplar differs';
  case 10, n ='Computer or nearest exemplar match category';
  otherwise, n='';
end
