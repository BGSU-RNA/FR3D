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
% Test with:  [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,BETA,UNITID,ModelNum,Readable,UseFile] = zReadCIFAtoms2('4V99.cifatoms')

% Filename may contain a path to the file; stripping away the path and extension leaves a potential PDB identifier

function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P,BETA,UNITID,ModelNum,Readable,UseFile] = zReadCIFAtoms(Filename,Verbose)

Verbose = 1;

CIFDownloaded = 0;
UseFile = '';

[pathstr, name, extension] = fileparts(Filename);
PDBID = name;

if exist(Filename),                    % file is found, either with a direct path, or on Matlab's path
  UseFile = Filename;
elseif exist([PDBID '.cif']) && length(pathstr) == 0,     % .cif file found, not directed to a specific path
  CIFDownloaded = 1;
  UseFile = [PDBID '.cif'];
elseif length(pathstr) > 0,            % directed to a specific path, but no file found there; should download cif to there

else                                   % no .cif file found
  if ~(exist('PDBFiles') == 7),        % if directory doesn't yet exist
    mkdir('PDBFiles');
  end
  path(path,[pwd filesep 'PDBFiles']);

  if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
    mkdir('PrecomputedData');
  end
  path(path,[pwd filesep 'PrecomputedData']);
end

if isempty(UseFile),
  try
    if Verbose > 0,
      fprintf('zReadCIFAtoms: Attempting to download %s.cif from PDB\n', PDBID);
      drawnow
    end
    if exist('webread'),
      c = webread(['http://www.rcsb.org/pdb/files/' PDBID '.cif']);
    else
      c = urlread(['http://www.rcsb.org/pdb/files/' PDBID '.cif']);
    end
    if Verbose > 0,
      fprintf('zReadCIFAtoms: Read the URL\n', PDBID);
      drawnow
    end
    fid = fopen(['PDBFiles' filesep PDBID '.cif'],'w');
    fprintf(fid,'%s\n',c);
    fclose(fid);
    CIFDownloaded = 1;
  catch
    fprintf('zReadCIFAtoms: Unable to download %s.cif from PDB\n', PDBID);
  end
end

if isempty(UseFile) && CIFDownloaded > 0 && ~isempty(strfind(lower(Filename),'.cifatoms')),
  try
    status = 1;
    GetFR3DPythonLocation
    if ~isempty(PythonLocation),
      fprintf('zReadCIFAtoms: Attempting to add unit ids and any crystal symmetries and save in %s\n',[PDBID '.cifatoms']);
      [status,result] = system([PythonVersion ' ' PythonLocation 'cifatom_writing.py ' pwd filesep 'PDBFiles' filesep PDBID '.cif']);
    else
      fprintf('zReadCIFAtoms: You can modify GetFR3DPythonLocation.m to specify where to find cifatom_writing.py to produce .cifatoms file.\n');
      fprintf('zReadCIFAtoms:  Python 2.7 is needed, also set the variable PythonVersion to tell how to run python 2.7 on your machine.\n');
      fprintf('repositories fr3d-python and pdbx are also needed.  Contact Craig Zirbel\n');
    end
    if status == 0,
      UseFile = [PDBID '.cifatoms'];
      fprintf('zReadCIFAtoms: Produced a .cifatoms file\n');
    else
      UseFile = [PDBID '.cif'];
      fprintf('zReadCIFAtoms: Was not able to produce a .cifatoms file\n');
      fprintf('%s\n',result)
    end
  catch
    fprintf('zReadCIFAtoms: Unable to use cifatom_writing.py to convert .cif file to .cifatoms file\n');
    UseFile = [PDBID '.cif'];
  end
end

DataSize = 1000;                       % keep track of how much space is allocated, allocate more as needed

if ~isempty(UseFile) && exist(UseFile),
  c = 1;                               % line counter
  if Verbose > 0,
    fprintf('zReadCIFAtoms: Reading %s\n', UseFile);
  end
  fid = fopen(UseFile,'r');
  StartHeader = 0;                       % keep track of whether any header lines have been found yet
  header = 1;
  headerlines = 0;
  fieldcounter = 0;
  row = 1;
  UnitIDField = [];
  ColumnsToKeep = 1;

  Line = fgetl(fid);
  while ischar(Line) && header == 1,
    headerlines = headerlines + 1;
    if StartHeader > 0 && any(strfind(Line,'ATOM') > 0) || any(strfind(Line,'HETATM') > 0),   % line contains ATOM or HETATM
      header = 0;                                    % header is over
    end
    if ~isempty(strfind(Line,'_atom_site.')),
      StartHeader = 1;
      fieldcounter = fieldcounter + 1;
      switch strrep(Line,' ',''),
      case '_atom_site.group_PDB'
        AtomTypeField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.id'
        AtomNumberField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.type_symbol'
      case '_atom_site.label_atom_id'
        AtomLabelIDField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_id'
      case '_atom_site.label_alt_id'
        VersionField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.label_comp_id'
        UnitNameField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.label_asym_id'
      case '_atom_site.label_entity_id'
      case '_atom_site.label_seq_id'
      case '_atom_site.pdbx_PDB_ins_code'
        InsertionCodeField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.Cartn_x'
        XCoordField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.Cartn_y'
        YCoordField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.Cartn_z'
        ZCoordField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
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
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.auth_comp_id'
      case '_atom_site.auth_asym_id'
        ChainField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.auth'
      case '_atom_id'
      case '_atom_site.pdbx_PDB_model_num'
        ModelNumberField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      case '_atom_site.unit_id'
        UnitIDField = fieldcounter;
        fieldtocolumn(fieldcounter) = ColumnsToKeep;
        ColumnsToKeep = ColumnsToKeep + 1;
      end
    end
    Line = fgetl(fid);
  end

  fclose(fid);

  headerlines = headerlines - 1;

  fprintf('zReadCIFAtoms: Found %d fields to read\n',fieldcounter);
  fprintf('zReadCIFAtoms: Found %d header lines to skip\n',headerlines);

  format = '';
  for i = 1:fieldcounter,
    if any(i == [XCoordField YCoordField ZCoordField ModelNumberField]),
      format = [format '%f'];
    elseif any(i == [AtomTypeField AtomNumberField AtomLabelIDField VersionField UnitNameField ChainField InsertionCodeField ResidueNumberField UnitIDField])
      format = [format '%s'];
    else
      format = [format '%*s'];            % skip this field
    end
  end

  if Verbose > 0,
    D = whos();
    fprintf('Bytes before reading:  %12d\n',sum(cat(1,D.bytes)));
  end

  fid = fopen(UseFile);
  A = textscan(fid,format,'headerlines',headerlines);
  fclose(fid);

  if Verbose > 0,
    D = whos();
    fprintf('Bytes  after reading:  %12d\n',round(sum(cat(1,D.bytes))));
  end

  i = 1;

  while strcmp(A{fieldtocolumn(AtomTypeField)}{i},'ATOM') || strcmp(A{fieldtocolumn(AtomTypeField)}{i},'HETATM'),
    i = i + 1;
  end
  NumLines = i-1;

  if Verbose > 0,
    D = whos();
    fprintf('%12d\n',round(sum(cat(1,D.bytes))));
  end

  fprintf('zReadCIFAtoms: Found %d lines of ATOM data\n',NumLines);

  ATOM_TYPE     = A{fieldtocolumn(AtomTypeField)}(1:NumLines);
  ATOMNAME      = A{fieldtocolumn(AtomLabelIDField)}(1:NumLines);
  ATOMNUMBER    = A{fieldtocolumn(AtomNumberField)}(1:NumLines);
  VERSION       = A{fieldtocolumn(VersionField)}(1:NumLines);
  UNITNAME      = A{fieldtocolumn(UnitNameField)}(1:NumLines);
  CHAIN         = A{fieldtocolumn(ChainField)}(1:NumLines);
  InsertionCode = A{fieldtocolumn(InsertionCodeField)}(1:NumLines);
  NTNUMBER      = A{fieldtocolumn(ResidueNumberField)}(1:NumLines);
  P             = [A{fieldtocolumn(XCoordField)}(1:NumLines) A{fieldtocolumn(YCoordField)}(1:NumLines) A{fieldtocolumn(ZCoordField)}(1:NumLines)];
  ModelNum      = A{fieldtocolumn(ModelNumberField)}(1:NumLines);

  if Verbose > 0,
    D = whos();
    fprintf('Bytes after storing: %12d\n',sum(cat(1,D.bytes)));
  end

  fprintf('zReadCIFAtoms: Defined most variables to pass back\n');

  BETA = NaN * ones(NumLines,1);
  OCC = cell(NumLines,1);
  for i = 1:NumLines,
    if InsertionCode{i} ~= '?'
      NTNUMBER{i} = [NTNUMBER{i} InsertionCode{i}];
    end
    ATOMNAME{i} = strrep(ATOMNAME{i},'"','');             % remove double quotes on backbone atom names
    OCC{i} = '';
  end

  if ~isempty(UnitIDField),
    UNITID = A{fieldtocolumn(UnitIDField)}(1:NumLines);
  else
    UNITID = OCC;
  end

  Readable = 1;

  if Verbose > 0,
    D = whos();
    fprintf('Bytes at the end: %12d\n',sum(cat(1,D.bytes)));
  end

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

