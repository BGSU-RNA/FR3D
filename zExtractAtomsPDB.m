% zExtractAtomsPDB(Filename) reads Filename.pdb and writes out lines 
% containing ATOM to Filename_Atoms.pdb

function [Header] = zExtractAtomsPDB(Filename,Outputname)

fid = fopen([Filename '.pdb'],'r');

Header.ModelStart = 1;

if fid > 0

  out = fopen([Outputname '.pdb'],'w');

  L = 1;

  c = 1;

  while L > -1
    L = fgets(fid);
    if L > -1
      if ~isempty(strfind(L,'RESOLUTION')),
        Header.Resolution = L;
      end
      if ~isempty(strfind(L,'EXPDTA')),
        Header.Expdata = L;
      end
      if strcmp(L(1:min(4,length(L))),'ATOM'),
        fprintf(out,'%s',L);
        c = c + 1;
        L = -1;                                   % stop reading header
      end
    end
  end

  L = 1;

  while L > -1
    L = fgets(fid);
    if L > -1
      if strcmp(L(1:min(4,length(L))),'ATOM'),
        fprintf(out,'%s',L);
        c = c + 1;
      end
      if strcmp(L(1:min(5,length(L))),'MODEL'),
        Header.ModelStart = [Header.ModelStart c];
      end
    end
  end

  fclose(fid);
  fclose(out);

  fprintf('Read %s.pdb for header information\n', Filename)

else

  fprintf('Could not open file %s.pdb\n', Filename);

end

if ~isfield(Header,'ExpData'),
  Header.ExpData = '';
end

if ~isfield(Header,'Resolution'),
  Header.Resolution = [];
end
