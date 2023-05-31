% zExemplarIDICalculationAll retrieves every exemplar combination and calculates the IDI between it and every other pair
% in the whole collection, writing the IDIs to the screen or a file

% This program could be written much more efficiently, but since it only needs to be run occasionally, it's OK for now.

Letters = 'ACGU';

if ~exist('Exemplar'),
  load PairExemplars
end

clear ExemplarIDI

fid = fopen('ExemplarIDIAll.txt','w');

for Class = [1:12],
  for Code1 = 1:4,
    for Code2 = 1:4,

      %fprintf('Calculating IDI values for %c%c %s\n', Letters(Code1), Letters(Code2), zEdgeText(Class,Code1,Code2));

      [NT1,NT2] = zGetExemplar(Class,Code1,Code2);

      if isempty(NT1.Code) || isempty(NT2.Code),
        IDI = NaN * ones(4,4);
      else
        NT1.ID = strrep(NT1.ID,'||.','');
        NT2.ID = strrep(NT2.ID,'||.','');

        if length(NT1.ID) == 0
          NT1.ID = 'NoID';
        end

        if length(NT2.ID) == 0
          NT2.ID = 'NoID';
        end

        for Class2 = [1:12]
          for a = 1:4,
            for b = 1:4,
              [M1,M2] = zGetExemplar(Class2,a,b);
              if isempty(M1.Code) || isempty(M2.Code),
                IDI(a,b) = NaN;
              else
                M1.ID = strrep(M1.ID,'||.','');
                M2.ID = strrep(M2.ID,'||.','');

                if length(M1.ID) == 0
                  M1.ID = 'NoID';
                end

                if length(M2.ID) == 0
                  M2.ID = 'NoID';
                end

                IDI(a,b) = zIsoDiscrepancy(NT1,NT2,M1,M2);

                fprintf('%s %c %c %s %s IDI %8.4f with %s %c %c %s %s ', zEdgeText(Class,Code1,Code2), Letters(Code1), Letters(Code2), NT1.ID, NT2.ID, IDI(a,b), zEdgeText(Class2,a,b), Letters(a), Letters(b), M1.ID, M2.ID);
                fprintf('http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s ', NT1.ID, NT2.ID);
                fprintf('http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s\n', M1.ID, M2.ID);

                fprintf(fid,'%s %c %c %s %s IDI %8.4f with %s %c %c %s %s ', zEdgeText(Class,Code1,Code2), Letters(Code1), Letters(Code2), NT1.ID, NT2.ID, IDI(a,b), zEdgeText(Class2,a,b), Letters(a), Letters(b), M1.ID, M2.ID);
                fprintf(fid,'http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s ', NT1.ID, NT2.ID);
                fprintf(fid,'http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s\n', M1.ID, M2.ID);

                IDI(a,b) = zIsoDiscrepancy(NT1,NT2,M2,M1);

                fprintf('%s %c %c %s %s IDI %8.4f with %s %c %c %s %s ', zEdgeText(Class,Code1,Code2), Letters(Code1), Letters(Code2), NT1.ID, NT2.ID, IDI(a,b), zEdgeText(-Class2,b,a), Letters(b), Letters(a), M2.ID, M1.ID);
                fprintf('http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s ', NT1.ID, NT2.ID);
                fprintf('http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s\n', M2.ID, M1.ID);

                fprintf(fid,'%s %c %c %s %s IDI %8.4f with %s %c %c %s %s ', zEdgeText(Class,Code1,Code2), Letters(Code1), Letters(Code2), NT1.ID, NT2.ID, IDI(a,b), zEdgeText(-Class2,b,a), Letters(b), Letters(a), M2.ID, M1.ID);
                fprintf(fid,'http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s ', NT1.ID, NT2.ID);
                fprintf(fid,'http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s\n', M2.ID, M1.ID);

              end
            end
          end
        end
      end

      ExemplarIDI{Class,Code1,Code2} = IDI;
      ExemplarIDI{100-Class,Code2,Code1} = IDI';

      %IDI

      %fprintf('Calculating IDI values for %c%c %s\n', Letters(Code2), Letters(Code1), zEdgeText(-Class,Code1,Code2));

      %IDI'

    end
  end
end

% save([pwd filesep 'FR3DSource' filesep 'PairExemplars'],'Exemplar','ExemplarIDI'); % Matlab version 7 only
% save PairExemplars_Version_6.mat Exemplar ExemplarIDI -V6 % for compatibility with older versions
