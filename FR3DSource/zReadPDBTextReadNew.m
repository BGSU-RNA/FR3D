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
    fprintf('zReadPDBTextReadNew: File does not consistently have 80 columns\n');
    Tcell = textread(Filename, '%s', 'delimiter', '\n');

    T = char(zeros(length(Tcell),80));     % character array of blanks
    keep = ones(length(Tcell),1);          % rows to keep, initially keep all

    for i = 1:length(Tcell),
      t = Tcell{i};
      a = min(length(t),80);
      T(i,1:a) = t(1:a);
      if a < 80,
%        keep(i) = 0;
      end
    end
    T = T(keep > 0,:);
  end
else
  if Verbose > 0,
    fprintf('zReadPDBTextReadNew: Attempting to download %s from PDB\n', Filename);
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
    T = [];
  end
end

if size(T,1) == 0,
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
%ATOM  A0001  N1    G 53064     241.043 256.642 248.664  1.00 79.27      25S  N

ATOM_TYPE  = strtrim(cellstr(T(k,1:6)));

ATOMNUMBER = str2num(T(k,7:11));
if isempty(ATOMNUMBER),
  ATOMNUMBER = 1:length(k);                 % not what the authors may have put, but fast
end

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
