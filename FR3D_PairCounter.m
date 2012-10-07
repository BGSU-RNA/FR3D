% zPairCountert loops through the molecules in File and counts how
% many instances of each pair type and each category occur, displaying the
% output in a table

% Load data like this:  File = zAddNTData('Nonredundant_list');

Classes = [];
PairCodes = [];

BClasses = [];
BPairCodes = [];

SClasses = [];
SPairCodes = [];

for f=1:length(File),
  for i = 1:length(File(f).Pair),
    if abs(File(f).Pair(i).Class) < 25,
      Classes = [Classes; [sprintf('%6.1f',File(f).Pair(i).Class) File(f).Pair(i).EdgeText ]];
      PairCodes = [PairCodes; File(f).Pair(i).Paircode];
    end
    if abs(File(f).Pair(i).Class) < 15,
      BClasses = [BClasses; [sprintf('%6.1f',File(f).Pair(i).Class) File(f).Pair(i).EdgeText ]];
      BPairCodes = [BPairCodes; File(f).Pair(i).Paircode];
    end
    if (abs(File(f).Pair(i).Class) > 15) && (abs(File(f).Pair(i).Class) < 25),
      SClasses = [SClasses; [sprintf('%6.1f',File(f).Pair(i).Class) File(f).Pair(i).EdgeText ]];
      SPairCodes = [SPairCodes; File(f).Pair(i).Paircode];
    end
  end
end

Labels = [];
[Table, Chi2, P, Labels] = crosstab(Classes,PairCodes);
for i=1:length(Labels(:,2)),
 if ~isempty(Labels{i,2}),
  switch str2num(Labels{i,2}),
    case 1, Labels{i,2} = 'AA';
    case 5, Labels{i,2} = 'AC';
    case 6, Labels{i,2} = 'CC';
    case 7, Labels{i,2} = 'GC';
    case 9, Labels{i,2} = 'AG';
    case 11, Labels{i,2} = 'GG';
    case 13, Labels{i,2} = 'AU';
    case 14, Labels{i,2} = 'CU';
    case 15, Labels{i,2} = 'GU';
    case 16, Labels{i,2} = 'UU';
  end
 end
end
zShowTable(Labels(:,1), Labels(:,2), Table, 'Basepairing and stacking');

Labels = [];
[Table, Chi2, P, Labels] = crosstab(BClasses,BPairCodes);
for i=1:length(Labels(:,2)),
 if ~isempty(Labels{i,2}),
  switch str2num(Labels{i,2}),
    case 1, Labels{i,2} = 'AA';
    case 5, Labels{i,2} = 'AC';
    case 6, Labels{i,2} = 'CC';
    case 7, Labels{i,2} = 'GC';
    case 9, Labels{i,2} = 'AG';
    case 11, Labels{i,2} = 'GG';
    case 13, Labels{i,2} = 'AU';
    case 14, Labels{i,2} = 'CU';
    case 15, Labels{i,2} = 'GU';
    case 16, Labels{i,2} = 'UU';
  end
 end
end
zShowTable(Labels(:,1), Labels(:,2), Table, 'Basepairing only');

Labels = [];
[Table, Chi2, P, Labels] = crosstab(SClasses,SPairCodes);
for i=1:length(Labels(:,2)),
 if ~isempty(Labels{i,2}),
  switch str2num(Labels{i,2}),
    case 1, Labels{i,2} = 'AA';
    case 5, Labels{i,2} = 'AC';
    case 6, Labels{i,2} = 'CC';
    case 7, Labels{i,2} = 'GC';
    case 9, Labels{i,2} = 'AG';
    case 11, Labels{i,2} = 'GG';
    case 13, Labels{i,2} = 'AU';
    case 14, Labels{i,2} = 'CU';
    case 15, Labels{i,2} = 'GU';
    case 16, Labels{i,2} = 'UU';
  end
 end
end
zShowTable(Labels(:,1), Labels(:,2), Table, 'Stacking only');
