% zReadCIFAtoms.m reads ATOM and HETATM records extracted from a CIF file
% It assumes that the data fields are always in this order:
% loop_
% _atom_site.group_PDB           'ATOM'
% _atom_site.id                  '699'
% _atom_site.type_symbol         'C'
% _atom_site.label_atom_id       'CD1'
% _atom_site.label_alt_id        '.'
% _atom_site.label_comp_id       'LEU'
% _atom_site.label_asym_id       '0'
% _atom_site.label_entity_id     '?'
% _atom_site.label_seq_id        '38'
% _atom_site.pdbx_PDB_ins_code   '?'
% _atom_site.Cartn_x             '-18.272'
% _atom_site.Cartn_y             '175.872'
% _atom_site.Cartn_z             '-187.65'
% _atom_site.occupancy           '?'
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
% _atom_site.auth_atom_id 
% _atom_site.pdbx_PDB_model_num 
% _atom_site.
%    'ATOM'    '699'    'C'    'CD1'    '.'    'LEU'    '0'    '?'    '38'    '?'    '-18.272'    '175.872'    '-187.65'    '?'    '?'    '?'    '?'    '?'    '?'
%    '?'    '.'    '38'    'LEU'    '0'    'CD1'    '1'    '2QBG|1|0|LEU|38'     ''

%   'ATOM'    '2'    'C'    'C'    '.'    'GLY'    '1'    '?'    '3'    '?'    '69.193'    '193.537'    '-102.468'    '?'    '?'    '?'    '?'    '?'    '?'    '?'
%    '.'    '3'    'GLY'    '1'    'C'    '1'    '2QBG|1|1|GLY|3'     ''

% Test with:  File = zAddNTData('2QBG.cifatoms')

function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,BETA,UNITID,ModelNum,Readable,UseFile] = zReadCIFAtoms(Filename,Verbose)

Verbose = 1;
CIFDownloaded = 0;
UseFile = '';

PDBID = strrep(upper(Filename),'.MAT','');
PDBID = strrep(PDBID,'.CIFATOMS','');
PDBID = strrep(PDBID,'.CIF','');
PDBID = strrep(PDBID,'-CIF','');

if exist(Filename),
  UseFile = Filename;
elseif exist([PDBID '.cif']),                           % no .cif file found
  CIFDownloaded = 1;
else
  try
    if Verbose > 0,
      fprintf('zReadCIFAtoms: Attempting to download %s.cif from PDB\n', PDBID);
    end
    c = urlread(['http://www.rcsb.org/pdb/files/' PDBID '.cif']);
    fid = fopen(['PDBFiles' filesep PDBID '.cif'],'w');
    fprintf(fid,'%s\n',c);
    fclose(fid);
    CIFDownloaded = 1;
  catch
    fprintf('zReadCIFAtoms: Unable to download %s.cif from PDB\n', PDBID);
  end
end

if CIFDownloaded > 0 && ~isempty(strfind(lower(Filename),'.cifatoms')),
  try
    if ~isempty(strfind(pwd,'zirbel')),
      fprintf('zReadCIFAtoms: Attempting to add unit ids and any crystal symmetries and save in .cifatoms file\n');
%      ['python C:\Users\zirbel\Documents\GitHub\fr3d-python\examples\cifatom_writing.py ' pwd filesep 'PDBFiles' filesep PDBID '.cif']
      system(['python C:\Users\zirbel\Documents\GitHub\fr3d-python\examples\cifatom_writing.py ' pwd filesep 'PDBFiles' filesep PDBID '.cif']);
    end
    UseFile = [PDBID '.cifatoms'];
  catch
    fprintf('zReadCIFAtoms: Unable to use cifatom_writing.py to convert .cif file to .cifatoms file\n');
  end
end

if ~isempty(UseFile),
  c = 1;                               % line counter
  if Verbose > 0,
    fprintf('zReadCIFAtoms: Reading %s\n', Filename);
  end

  T = textread(Filename,'%s','whitespace','\n');

  fprintf('zReadCIFAtoms: Read %d lines from file\n',size(T,1));

  header = 1;
  fieldcounter = 0;  
  row = 1;
  UnitIDField = [];
  while header == 1,
    if any(strfind(T{row,1},'ATOM') == 1) || any(strfind(T{row,1},'HETATM') == 1),                % line starts with ATOM or HETATM
      header = 0;                                    % header is over
      break
    end
    if ~isempty(strfind(T{row,1},'_atom_site.')),
      fieldcounter = fieldcounter + 1;
      switch strrep(T{row,1},' ',''),
      case '_atom_site.group_PDB' 
        AtomTypeField = fieldcounter;
      case '_atom_site.id'
        AtomNumberField = fieldcounter;
      case '_atom_site.type_symbol'
      case '_atom_site.label_atom_id'
        AtomLabelIDField = fieldcounter;
      case '_atom_id' 
      case '_atom_site.label_alt_id' 
        VersionField = fieldcounter;
      case '_atom_site.label_comp_id' 
        UnitNameField = fieldcounter;
      case '_atom_site.label_asym_id' 
        ChainField = fieldcounter;
      case '_atom_site.label_entity_id' 
      case '_atom_site.label_seq_id' 
      case '_atom_site.pdbx_PDB_ins_code' 
        InsertionCode = fieldcounter;
      case '_atom_site.Cartn_x' 
        XCoordField = fieldcounter;
      case '_atom_site.Cartn_y' 
        YCoordField = fieldcounter;
      case '_atom_site.Cartn_z' 
        ZCoordField = fieldcounter;
      case '_atom_site.occupancy' 
      case '_atom_site.B_iso_or_equiv' 
      case '_atom_site.Cartn_x_esd' 
      case '_atom_site.Cartn_y_esd' 
      case '_atom_site.Cartn_z_esd' 
      case '_atom_site.occupancy_esd' 
      case '_atom_site.B_iso_or_equiv_esd' 
      case '_atom_site.pdbx_formal_charge' 
      case '_atom_site.auth_seq_id'
        ResidueNumberField = fieldcounter;
      case '_atom_site.auth_comp_id' 
      case '_atom_site.auth_asym_id' 
      case '_atom_site.auth'
      case '_atom_id' 
      case '_atom_site.pdbx_PDB_model_num'
        ModelNumberField = fieldcounter;
      case '_atom_site.unit_id'
        UnitIDField = fieldcounter;
      end 
    end
    row = row + 1;
  end

  AtomCounter = 1;

  while row <= size(T,1),
    entry = zStringSplit(T{row,1});

    if strcmp(entry{1},'ATOM') || strcmp(entry{1},'HETATM'),

%function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,OCC,BETA,ModelNum,Readable] = zReadCIFAtoms(Filename,Verbose)

      ATOM_TYPE{AtomCounter} = entry{AtomTypeField};
      ATOMNAME{AtomCounter} = strrep(entry{AtomLabelIDField},'"','');       % remove double quotes on backbone atom names
      ATOMNUMBER{AtomCounter} = entry{AtomNumberField};
      VERSION{AtomCounter} = entry{VersionField};
      UNITNAME{AtomCounter} = entry{UnitNameField};
      CHAIN{AtomCounter} = entry{ChainField};
      if entry{InsertionCode} ~= '?'
        NTNUMBER{AtomCounter} = [entry{ResidueNumberField} entry{InsertionCode}];
      else
        NTNUMBER{AtomCounter} = entry{ResidueNumberField};
      end
      P(AtomCounter,:) = [str2num(entry{XCoordField}) str2num(entry{YCoordField}) str2num(entry{ZCoordField})];
      BETA(AtomCounter,1) = NaN;
      OCC{AtomCounter} = 'NaN';
      ModelNum(AtomCounter,1) = str2num(entry{ModelNumberField});
      if isempty(UnitIDField),
        UNITID{AtomCounter} = '';
        MN = 1;                          % just in case
      else
        UNITID{AtomCounter} = entry{UnitIDField};
        UnitIDParts = zStringSplit(entry{UnitIDField},'|');
        MN = str2num(UnitIDParts{2});
        if MN ~= ModelNum(AtomCounter,1),
          fprintf('zReadCIFAtoms: Disagreement about model number in this entry: %s\n',T{row,1});
        end
      end

    end
    AtomCounter = AtomCounter + 1;
    row = row + 1;

  end

  Readable = 1;

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
    UNITID   = [];
    Readable = 0;
    UseFile = '';
end

return









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

