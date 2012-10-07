% zWriteExemplarPDB reads the file of pair exemplars and writes them out
% to a single PDB file, spaced 20 Angstroms apart in a plane

function [void] = zWriteExemplarPDB(g)

% load exemplars -------------------------------------

load('PairExemplars','Exemplar');

if nargin < 1,

  % loop through pairs and classifications ----------------------

  fid = fopen('Isostericity\PairExemplarPDB.pdb','w');       % open for writing

  a = 1;                                         % atom number
  c = 1;                                         % pair number

  for pc = 1:length(Exemplar(1,:)),
    for row = 1:length(Exemplar(:,pc)),

      E = Exemplar(row,pc);
  
      if ~isempty(E.NT1),
        R = E.NT1.Rot;
        sh = E.NT1.Center;

        a = zWriteNucleotidePDB(fid,E.NT1,a,c,R,sh);
        a = zWriteNucleotidePDB(fid,E.NT2,a,c,R,sh);
        c = c + 1;
      end

    end
  end

  fclose(fid);

  fprintf('Wrote PairExemplarPDB.pdb\n');

else

  for pc = 1:length(Exemplar(1,:)),
    for row = 1:length(Exemplar(:,pc)),

      a = 1;                                         % atom number

      E = Exemplar(row,pc);
  
      if ~isempty(E.NT1),
        R  = E.NT1.Rot';
        sh = E.NT1.Fit(1,:);

        s = abs(E.Class) - fix(abs(E.Class));
        if s > 0,
          S = ['_Subcat' num2str(10*s)];
        else
          S = '';
        end

        fid = fopen(['Isostericity\PairExemplar_' strrep(E.Pair.EdgeText,' ','') '_' E.NT1.Base E.NT2.Base '_' E.Filename '_' E.NT1.Number E.NT1.Chain '_' E.NT2.Number E.NT2.Chain S '.pdb'],'w');       % open for writing

        a = zWriteNucleotidePDB(fid,E.NT1,a,0,R,sh);
        a = zWriteNucleotidePDB(fid,E.NT2,a,0,R,sh);

        fclose(fid);
      end

    end
  end


end
