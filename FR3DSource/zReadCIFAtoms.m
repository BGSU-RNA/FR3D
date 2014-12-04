% zReadCIFAtoms.m reads ATOM and HETATM records extracted from a CIF file
% It assumes that the data fields are always in this order:
% loop_
% _atom_site.group_PDB 
% _atom_site.id 
% _atom_site.type_symbol 
% _atom_site.label% _atom_id 
% _atom_site.label_alt_id 
% _atom_site.label_comp_id 
% _atom_site.label_asym_id 
% _atom_site.label_entity_id 
% _atom_site.label_seq_id 
% _atom_site.pdbx_PDB_ins_code 
% _atom_site.Cartn_x 
% _atom_site.Cartn_y 
% _atom_site.Cartn_z 
% _atom_site.occupancy 
% _atom_site.B_iso_or_equiv 
% _atom_site.Cartn_x_esd 
% _atom_site.Cartn_y_esd 
% _atom_site.Cartn_z_esd 
% _atom_site.occupancy_esd 
% _atom_site.B_iso_or_equiv_esd 
% _atom_site.pdbx_formal_charge 
% _atom_site.auth_seq_id 
% _atom_site.auth_comp_id 
% _atom_site.auth_asym_id 
% _atom_site.auth% _atom_id 
% _atom_site.pdbx_PDB_model_num 

% Test with:  File = zAddNTData('2QBG_temp.cifatoms')

function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,OCC,BETA,ModelNum,Readable] = zReadCIFAtoms(Filename,Verbose)

Verbose = 1;

KeepPDB = 1;

if exist(Filename),
  c = 1;                               % line counter
  if Verbose > 0,
    fprintf('zReadCIFAtoms: Reading %s\n', Filename);
  end

  T = textread(Filename,'%s','whitespace','\n');

  fprintf('zReadCIFAtoms:  Read %d lines from file\n',length(Lines));

  keyboard

else
  fprintf('zReadCIFAtoms: File %s was not found\n', Filename);
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

