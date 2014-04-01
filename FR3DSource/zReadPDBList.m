% zReadPDBList(Filename) reads Filename.pdb for a list of PDB file names

function [NewNames] = zReadPDBList(Filename,Verbose)

if nargin < 2,
  Verbose = 1;
end

NewNames = {};

if strcmpi(Filename,'AllFiles_list'),
  [s,NewNames] = mGetPDBFilenames;
elseif ~isempty(strfind(Filename,'http://')),
  try
    c = urlread(Filename);
  catch
    fprintf('Not able to retrieve a list of filenames from the URL %s\n',Filename);
  end
  try 
    lines = zStringSplit(c,char(10));
    for i = 1:length(lines),
      m = zStringSplit(lines{i},',');
      if length(m) >= 3,
        fn = strrep(m{2},'"','');
        if ~isempty(fn),
          NewNames{i} = fn;
        end
      end
    end
  catch
    fprintf('Something wrong with the list of filenames from the URL %s\n',Filename);    
  end
  try
    fn = strrep(Filename,'http://rna.bgsu.edu/rna3dhub/nrlist/download/','');
    fn = strrep(fn,'/csv','');
    fn = strrep(fn,'/','_');
    fn = ['nrlist_' fn '_list.pdb'];
    if length(NewNames) > 0,
      fid = fopen(['PDBFiles' filesep fn],'w');
      for i = 1:length(NewNames),
        fprintf(fid,'%s\n',NewNames{i});
      end
      fclose(fid);
      fprintf('Wrote downloaded file list %s to the folder PDBFiles\n',fn);
    end
  end
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

    L = '';

    while ischar(L),
      L = fgetl(fid);
      if ischar(L),
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
