% zEnterPairSelection(P) is part of zPairViewer.  It prompts the user to
% specify the pair code, category, group, and/or box to display.  Current
% values are given by the structured variable P; this is returned when
% modified.

function [P] = zEnterPairSelection(P)

% --------- Set default values for selection parameters ---------------

if isempty(P),
  P.Paircode = 13; 
  P.Category  = 1;
  P.Decimal   = 0;        % 1 - use 1.0 only; 0 - round to 1
  P.Group     = 1;        % hand, computer, both, etc.
  P.Sequential= 0;
  P.Expert    = 0;
  P.Inbox     = [-5 5 -5 5 -5 5 -1.1 1.1 -95 275];  % almost everything
  P.Context   = 3;
  P.PairFilename = '1s72';
  P.PairNucleotide1 = '11';
  P.PairNucleotide2 = '12';
  P.PairDisc  = 5;
end

% --------- Prompt user to change paircode ---------------

fprintf('\n');
fprintf('1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU\n')
fprintf('Enter paircode(s) to view        [');
fprintf(' %2d',P.Paircode);
fprintf('] ');
inp = input('');
if ~isempty(inp),
  P.Paircode = inp;
end

% --------- Prompt user to change category ---------------

fprintf('\n');
if P.Expert == 0,
  fprintf('  0  Display all categories together\n');
  zListCategoryNames;
end

if P.Decimal == 1,
  fprintf('Enter category(ies) to view     [');
  fprintf('%1.2f ',P.Category);
  fprintf('] ');
else
  fprintf('Enter category(ies) to view     [');
  fprintf('%1d ',P.Category);
  fprintf('] ');
end

inp = lower(input('','s'));
if ~isempty(str2num(inp)),
  P.Category = str2num(inp);
  P.Decimal = ~isempty(strfind(inp,'.'));
end

if any(P.Category == 50),
  fprintf('\n');
  fprintf('Specify [xmin xmax ymin ymax zmin zmax n3min n3max angmin angmax]\n');
  fprintf('Current value is (hit return to use it)\n');
  fprintf('[ ');
  fprintf('%7.2f', P.Inbox);
  fprintf(' ]\n');
  inp = input('');
  if ~isempty(inp)
    P.Inbox = inp;
  end
end

if any(P.Category == 51),
  fprintf('\n');
  fprintf('Enter file name as text [%s] ',P.PairFilename);
  inp = input('','s');
  if ~isempty(inp),
    P.PairFilename = inp;
  end
  fprintf('Enter first nucleotide number  [%s] ',P.PairNucleotide1);
  inp = input('','s');
  if ~isempty(inp),
    P.PairNucleotide1 = inp;
  end
  fprintf('Enter second nucleotide number [%s] ',P.PairNucleotide2);
  inp = input('','s');
  if ~isempty(inp),
    P.PairNucleotide2 = inp;
  end
  fprintf('Enter maximum acceptable discrepancy [%d] ',P.PairDisc);
  inp = input('');
  if ~isempty(inp),
    P.PairDisc = inp;
  end
end

% --------------- Prompt user to change group -----------------------------

fprintf('\n');
if ~any(P.Category == 0) & all(P.Category < 40),
  if P.Expert == 0,
    for k=1:10,
      fprintf('Group %1d - %s\n',k,zGroupNames(k));
    end  
  end

  fprintf('Enter group to view        [%1d] ',P.Group);
  inp = input('');
  if ~isempty(inp),
    P.Group = inp(1);
  end
end

if any(P.Category == 0) | (P.Group == 6),
  fprintf('Enter 0 for all pairs, 1 for sequential separated by 1, etc. ');
  fprintf('[%1d] ', P.Sequential);
  inp = input('');
  if ~isempty(inp),
    if inp == 0,
      P.Group = 1;                           % doesn't matter
      P.Sequential = 0;
    else
      P.Group = 6;
      P.Sequential = inp;
    end
  else
    if P.Sequential > 0,
      P.Group = 6;
    end
  end
else
  P.Sequential = 0;
end    
