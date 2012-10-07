% zReadPDBList(Filename) reads Filename.pdb for a list of PDB file names

function [NewNames] = zReadPDBList(Filename)

NewNames = '';

if isempty(strfind(Filename,'_list')),
  NewNames = {Filename};
elseif strcmp(Filename,'AllFiles_list'),
  [s,NewNames] = mGetPDBFilenames;
elseif ~isempty(Filename),
  fid = fopen(['PDBFiles' filesep Filename '.pdb'],'r');

  if fid > 0

    L = 1;

    while L > -1
      L = fgetl(fid);
      if L > -1
        if ~isempty(strfind(L,'_list')),
          NewNames = [NewNames; zReadPDBList(L)];
        else
          NewNames = [NewNames; {L}];
        end
      end
    end

    fclose(fid);

    fprintf('Read %s.pdb\n', Filename)

  else

    fprintf('Could not open file %s.pdb\n', Filename);

  end

end
