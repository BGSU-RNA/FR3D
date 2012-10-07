% some PDB files have some atoms with two "versions" A and B.

function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,OCC,BETA,ModelNum,Readable] = zReadPDBTextReadNew(Filename,Verbose)

KeepPDB = 1;

if exist(Filename),
  if Verbose > 0,
    fprintf('Reading %s\n', Filename);
  end
  try
    T = textread(Filename,'%80c');
  catch
    T = [];
    fid = fopen(Filename,'r');             % make all lines 80 characters
    if fid > 0
      L = 1;
      while L > -1
        L = fgets(fid);
        if L > -1
          a = length(L);
          if a < 80,
            L = [L ' '*ones(1,80-a)];
          elseif a > 80,
            L = L(1:80);
          end
          T = [T; L];
        end
      end
    end
  end
else
  if Verbose > 0,
    fprintf('Attempting to download %s from PDB\n', Filename);
  end
  a = ['http://www.rcsb.org/pdb/files/' Filename];
  try
    c = urlread(a);
    L = length(c);
    T = reshape(c',81,[])';
    T = T(:,1:80);

    fid = fopen(['PDBFiles' filesep Filename],'w');
    for i = 1:length(T(:,1)),
      fprintf(fid,'%s\n',T(i,:));
    end
    fclose(fid);

  catch
    fprintf('Unable to download %s from PDB\n', Filename);
    ATOM_TYPE = [];
    ATOMNUMBER = [];
    NTLETTER = [];
    NTNUMBER = [];
    P        = [];
    ATOMNAME = [];
    VERSION  = [];
    UNITNAME = [];
    OCC      = [];
    BETA     = [];
    CHAIN    = [];
    ModelNum = [];
    Readable = 0;
    return
  end
end

% -------------------------------------- Extract MODRES rows

j = find( (T(:,1) == 'M') .* (T(:,2) == 'O') .* (T(:,3) == 'D') .* (T(:,4) == 'R') .* (T(:,5) == 'E') .* (T(:,6) == 'S'));

if 0 > 1,
  fid = fopen('ModifiedNucleotides.txt','a');
  for i = 1:length(j),
    if T(j(i),26) == ' ',                     % avoid modified am. acid or DNA
      fprintf(fid,'%s\n',T(j(i),:));
    end
  end
  fclose(fid);
end

%         1         2         3         4         5         6         7         
%1234567890123456789012345678901234567890123456789012345678901234567890123456789
%MODRES 1H3E PSU B   35    U  PSEUDOURIDINE-5'-MONOPHOSPHATE                  
%MODRES 1H3E 5MU B   54    T  RIBOSYLTHYMINE-5'-MONOPHOSPHATE                  
%MODRES 1H3E PSU B   55    U  PSEUDOURIDINE-5'-MONOPHOSPHATE                   
%MODRES 1H3E MAD B   58    A  6-HYDRO-1-METHYLADENOSINE-5'-MONOPHOSPHATE       

% -------------------------------------- Extract ATOM, HETATM, MODEL rows

i = find( (T(:,1) == 'A') .* (T(:,2) == 'T') .* (T(:,3) == 'O') .* (T(:,4) == 'M'));

j = find( (T(:,1) == 'H') .* (T(:,2) == 'E') .* (T(:,3) == 'T') .* (T(:,4) == 'A') .* (T(:,5) == 'T') .* (T(:,6) == 'M'));

k = [i; j];

k = sort(k);

m = find( (T(:,1) == 'M') .* (T(:,2) == 'O') .* (T(:,3) == 'D') .* (T(:,4) == 'E') .* (T(:,5) == 'L'));

if isempty(m),
  m = 0;
end

for z = 1:length(k),
  ModelNum(z,1) = length(find(k(z) > m));
end

%         1         2         3         4         5         6         7         
%1234567890123456789012345678901234567890123456789012345678901234567890123456789
%ATOM     20  OP2   A 0  11      20.863 146.760 101.535  1.00 65.62           O 
%HETATM13211  OP1 1MA 0 628      40.887 106.997  94.641  1.00 24.89           O 
%HETATM91248  O   HOH     1      82.803  66.117  95.494  1.00 32.08           O 
%ATOM   3047  C6    C A 141A     98.940 267.136 -49.615  1.00 38.03           C 

ATOM_TYPE  = strtrim(cellstr(T(k,1:6)));
ATOMNUMBER = str2num(T(k,7:11));
ATOMNAME   = strtrim(cellstr(T(k,12:16)));
VERSION    = cellstr(T(k,17));
UNITNAME   = strtrim(cellstr(T(k,18:20)));
CHAIN      = strtrim(cellstr(T(k,21:22)));  % some entries are blank, for water
NTNUMBER   = strtrim(cellstr(T(k,23:27)));  % allow for 141A, etc.
X          = str2num(T(k,31:38));
Y          = str2num(T(k,39:46));
Z          = str2num(T(k,47:54));

P          = [X Y Z];
OCC        = str2num(T(k,55:60));
BETA       = str2num(T(k,61:66));

Readable = 1;

