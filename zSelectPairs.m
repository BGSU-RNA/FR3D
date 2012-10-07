% zSelectPairs(File,Param) returns a list of indices of selected bases and
% pairs using selection parameters in Param.
% Its output is an array SP of structured variables.
% For instance, SP(1) contains these fields:
%   SP(1).Filenum        The number of the file from which this pair comes
%   SP(1).B1Index        The index of the first base in the pair
%   SP(1).B2Index        The index of the second base in the pair
%   SP(1).PairIndex      The index of the pair of bases
%   SP(1).CI             The hand index of the pair (may be zero)

function [SP] = zSelectPairs(File,Param)

k = 1;                                    % index of kept pairs

if any(Param.Category == 51),             % match a given pair
  Pf = zFileNumberLookup(File,Param.PairFilename);  % must already be loaded
  if Pf > 0,
    PIndex1 = zIndexLookup(File(Pf),{Param.PairNucleotide1});
    PIndex2 = zIndexLookup(File(Pf),{Param.PairNucleotide2});
    MatchPair =zAnalyzePair(File(Pf).NT(PIndex1(1)),File(Pf).NT(PIndex2(1)));
  end
end

for f = 1:length(File),                   % Loop through each file
 for i = 1:length(File(f).Pair)           % Loop through pairs

  keep = 0;                               % Default is not to keep this pair
  p = File(f).Pair(i);                    % Store this pair separately

  ci = File(f).CI(p.Base1Index,p.Base2Index);  % Comment index
  if ci > 0,                                   % Pair is in hand file
    HandClass = File(f).HandClass(ci);         % Use the hand class
  else
    HandClass = 0;                             % Hand class 0
  end

  % -------- Match paircode ----------------------------------------------

  if any(Param.Paircode == p.Paircode),      % if a paircode matches

  % -------- Round categories to nearest integer, if desired ------------- 

  if Param.Decimal == 1,                  % select with decimal places in categ
    pClass = p.Class;
    hClass = HandClass; 
    eClass = p.Classes(1);
  else
    pClass = fix(p.Class);                % remove decimal part
    hClass = fix(HandClass);
    eClass = fix(p.Classes(1));
  end

  PC = Param.Category;                    % desired category(ies)

  % -------- Restrict to indices with desired category and group

 if any(p.Base2Index - p.Base1Index == Param.Sequential) ...
     & (File(f).NT(p.Base1Index).Chain == File(f).NT(p.Base2Index).Chain)
   Sequential = 1;
 else
   Sequential = 0;
 end;

 switch Param.Group,                       
  case 1,
    if any(pClass == PC),           % match computer category only
      keep = 1;
    end
  case 2,
    if any(hClass == PC),           % match hand category only
      keep = 1;
    end
  case 3,
    if any(pClass == PC) | any(hClass == PC), % match either
      keep = 1;
    end
  case 4,
    if any(pClass == PC) & ~(hClass == pClass),  % computer, hand differs
      keep = 1;
    end
  case 5,
    if any(hClass == PC) & ~(hClass == pClass),  % hand, computer differs
      keep = 1;
    end
  case 6,
    if any(pClass == PC) & (Sequential == 1),
      keep = 1;
    end
  case 7,
    if any(eClass == PC) & (p.Distances(1) < 5),  % match nearest exemplar only
      keep = 1;
    end
  case 8,
    if any(eClass == PC) & ~(eClass == pClass) & (p.Distances(1) < 5),  % match exemplar, not cutoff
      keep = 1;
    end
  case 9,
    if any(pClass == PC) & ~(eClass == pClass) & (p.Distances(1) < 5), % match cutoff, not exemplar
      keep = 1;
    end
  case 10,
    if any(pClass == PC) | (any(eClass == PC) & (p.Distances(1) < 5)), % match cutoff or exemplar
      keep = 1;
    end
  case 11,
    if any(bClass == PC) & (p.Distances(1) < 5),  % match nearest exemplar only
      keep = 1;
    end
  case 12,
    if any(bClass == PC) & ~(bClass == pClass) & (p.Distances(1) < 5),  % match exemplar, not cutoff
      keep = 1;
    end
  case 13,
    if any(pClass == PC) & ~(bClass == pClass) & (p.Distances(1) < 5), % match cutoff, not exemplar
      keep = 1;
    end
  case 14,
    if any(pClass == PC) | (any(bClass == PC) & (p.Distances(1) < 5)), % match cutoff or exemplar
      keep = 1;
    end
 end

 for j=1:length(Param.Category)           % go through all desired categories
  switch Param.Category(j),
    case  0, if (Param.Group < 6) | ...
                ((Param.Group == 6) & ...
                 (Sequential == 1)),
               keep = 1;                            % display all pairs
             end               
    case 40, inbox = [-12 12 -12 12  2  5  0   1.1 -95 275]; % not very good!
    case 41, inbox = [-12 12 -12 12  2  5 -1.1 0   -95 275];
    case 42, inbox = [-12 12 -12 12 -5 -2  0   1.1 -95 275];
    case 43, inbox = [-12 12 -12 12 -5 -2 -1.1 0   -95 275];
    case 44, inbox = [-12 12 -12 12 -2  2 -1.1 1.1 -95 275];
    case 45, inbox = [-12 12 -12 12 -1  1 -1.1 1.1 -95 275];
    case 50, inbox = Param.Inbox;
  end

  if any(Param.Category(j) == [40 41 42 43 44 45 50]),     % check inbox
    if Param.Category(j) <= 45,
     n1 = File(f).NT(p.Base1Index);             % first base
     n2 = File(f).NT(p.Base2Index);             % second base
    end
    if (p.Displ(:,1)  > inbox(1)) ...
     & (p.Displ(:,1)  < inbox(2)) ...
     & (p.Displ(:,2)  > inbox(3)) ...
     & (p.Displ(:,2)  < inbox(4)) ...
     & (p.Displ(:,3)  > inbox(5)) ...
     & (p.Displ(:,3)  < inbox(6)) ...
     & (p.Normal(:,3) > inbox(7)) ...
     & (p.Normal(:,3) < inbox(8)) ...
     & (p.Ang(:)      > inbox(9)) ...
     & (p.Ang(:)      < inbox(10)) ...
     & ((Param.Category(j) > 43) | ((abs(p.Gap) > 2.0) & (p.MinDist < 3.5))) ...
     & ((Param.Category(j) ~= 44) | (abs(p.Gap) < 1.2)) ...
     & ((Param.Category(j) ~= 45) | ((abs(p.Gap) < 1.0) & (p.MinDist < 3) ...
                                  & (abs(p.Normal(3)) > 0.7) & (abs(p.Class) >= 30))), ...
     keep = 1;
    end
  end
 end


 % choice 51 is flawed:  it only looks at pairs that are already in File.Pair,
 % but those all have i < j, so if you want to match with i > j, it won't
 % happen.

 if any(Param.Category == 51),
   pd = zPairDiscrepancy(MatchPair,p);
   if pd < Param.PairDisc,
     keep = 1;
   end
 else
   pd = 0;
 end 

 if (keep == 1),
   n.Filenum   = f;            % number of the file from which this pair comes
   n.B1Index   = p.Base1Index;   % The index of the first base in the pair
   n.B2Index   = p.Base2Index;   % The index of the second base in the pair
   n.PairIndex = i;              % The index of the pair of bases
   n.CI        = ci;             % The hand index of the pair (may be zero)
   n.HandClass = HandClass;
   n1          = File(f).NT(p.Base1Index);             % first base
   n2          = File(f).NT(p.Base2Index);             % second base
   n.MinDist   = p.MinDist;
   n.C1pC1p    = zDistance(n1.Sugar(1,:),n2.Sugar(1,:));  % C1' - C1' distance
   n.PairDisc  = pd;
   SP(k) = n;                    % Append this pair to the selected pairs (SP)
   k = k + 1;
 end

end  % if paircode matches
end  % loop through i
end  % loop through f

Pair = {'AA'; 'CA'; 'GA'; 'UA'; 'AC'; 'CC'; 'GC'; 'UC'; 'AG'; 'CG'; ...
        'GG'; 'UG'; 'AU'; 'CU'; 'GU'; 'UU'};

%fprintf('\n');

if ~exist('SP'),
  fprintf('No ');
  for i=1:length(Param.Paircode),
    fprintf(' %s ', Pair{Param.Paircode(i)});
  end
  for i=1:length(Param.Category),
    fprintf(' %s ', zCategoryName(Param.Category(i)));
  end
  fprintf(' pairs were found.\n');
  SP = [];
else
  if Param.Category(1) == 0,
    fprintf('There are %1d ', length(SP));
    for i=1:length(Param.Paircode),
      fprintf(' %s ', Pair{Param.Paircode(i)});
    end
    fprintf(' pairs.\n');
  else
    fprintf('Found %6d ', length(SP));
    for i=1:length(Param.Paircode),
      fprintf(' %s ', Pair{Param.Paircode(i)});
    end
    for i=1:length(Param.Category),
      fprintf(' %s ', zCategoryName(Param.Category(i)));
    end
    fprintf(' pairs.\n');
  end          
end

 if any(Param.Category == 51),
   VP.SortKeys = 23;
   SP = zSortPairs(File,SP,VP);  % sort by distance to the given pair
 end 
