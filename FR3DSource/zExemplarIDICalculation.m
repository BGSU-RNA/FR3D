% zExemplarIDICalculation retrieves every exemplar combination and calculates the IDI between it and every other pair in the family, saving the IDIs in the variable ExemplarIDI in PairExemplars.mat

% This program could be written much more efficiently, but since it only needs to be run occasionally, it's OK for now.

Letters = 'ACGU';

if ~exist('Exemplar'),
  load PairExemplars
end

clear ExemplarIDI

for Class = [1:15],
  for Code1 = 1:4,
    for Code2 = 1:4,

fprintf('Calculating IDI values for %c%c %s\n', Letters(Code1), Letters(Code2), zEdgeText(Class,Code1,Code2));

[NT1,NT2] = zGetExemplar(Class,Code1,Code2);

if isempty(NT1.Code) || isempty(NT2.Code),
  IDI = NaN * ones(4,4);
else
  for a = 1:4,
    for b = 1:4,
      [M1,M2] = zGetExemplar(Class,a,b);
      if isempty(M1.Code) || isempty(M2.Code),
        IDI(a,b) = NaN;
      else
        IDI(a,b) = zIsoDiscrepancy(NT1,NT2,M1,M2);
      end
    end
  end
end

ExemplarIDI{Class,Code1,Code2} = IDI;
ExemplarIDI{100-Class,Code2,Code1} = IDI';

IDI

fprintf('Calculating IDI values for %c%c %s\n', Letters(Code2), Letters(Code1), zEdgeText(-Class,Code1,Code2));

IDI'

    end
  end
end

save([pwd filesep 'FR3DSource' filesep 'PairExemplars'],'Exemplar','ExemplarIDI'); % Matlab version 7 only
save PairExemplars_Version_6.mat Exemplar ExemplarIDI -V6 % for compatibility with older versions
