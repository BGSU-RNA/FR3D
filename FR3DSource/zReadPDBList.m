% zReadPDBList(Filename) reads Filename.pdb for a list of PDB file names

function [NewNames] = zReadPDBList(Filename,Verbose)

if nargin < 2,
  Verbose = 1;
end

NewNames = '';

if strcmpi(Filename,'AllFiles_list'),
  [s,NewNames] = mGetPDBFilenames;
elseif ~isempty(strfind(Filename,'_equiv')),
  load PDBInfo
  Reference = strrep(Filename,'_equiv','');
  NewNames  = {Reference};
  for i = 1:length(t(:,1)),                      % loop through PDB files
    if strcmpi(t{i,10},Reference),
      NewNames = [NewNames; t(i,1)];
    end
  end
elseif isempty(strfind(Filename,'_list')),
  NewNames = {Filename};
elseif ~isempty(Filename),
  Filename = strrep(Filename,'.pdb','');
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

    fprintf('Read list %s.pdb\n', Filename)

  else

    fprintf('Could not open file %s.pdb\n', Filename);

  end

end
