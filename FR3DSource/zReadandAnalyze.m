% zReadandAnalyze.m reads a .pdb or .cif file, then
% pulls out the locations of base atoms,
% and computes the best shift and rotation from a standard base

% It returns a variable File with fields:
%
% NT        array of nucleotide structures
% NumNT     the total number of distinct nucleotides found

% Each nucleotide has these fields:
%  Base     Nucleotide letter A, C, G, U
%  Code     Numeric code for nucleotide letter, A=1, C=2, G=3, U=4
%  Chain    Which chain of the molecule the base is located in
%  Number   Cell array of nucleotide number within the molecule
%  Loc      Observed locations of base atoms
%  Center   Geometric center of observed base atoms
%  Sugar    Observed locations of sugar atoms
%  Rot      3x3 matrix giving rotation from standard to observed locations
%  Scale    Best scaling from standard to observed locations
%  Shift    Best shift from standard to observed locations
%  Fit      Best fit of observed locations to standard base locations

% Filename may contain a path to the file; stripping away the path and extension leaves a potential PDB identifier

function [File] = zReadandAnalyze(Filename,Verbose)

[pathstr, PDBID, extension] = fileparts(Filename);

if nargin < 2,
  Verbose = 1;
end

% read Filename ----------------------------------------------------

Stop = 0;

if exist([pwd filesep 'PDBFiles' filesep 'Trouble reading'], 'dir') == 7
  if exist([pwd filesep 'PDBFiles' filesep 'Trouble reading' filesep PDBID], 'file') == 2
    Stop = 1;
    fprintf('zReadandAnalyze is skipping %s because it appeared in FR3D/PDBFiles/Trouble Reading\n', PDBID);
    return
  end
end

if ~isempty(strfind(lower(extension),'.pdb'))
  [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P, OCC, BETA, ModelNum, Readable] = zReadPDBTextReadNew(Filename,Verbose);
  CoordinateFile = Filename;
  if Verbose > 0
    fprintf('zReadandAnalyze: Read PDB file %s\n',Filename);
  end
else
  % read a special version of the .cif file that has symmetry operators applied and unit ids supplied
  [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, UNITNAME, CHAIN, NTNUMBER, P, BETA, UNITID, ModelNum, Readable, CoordinateFile] = zReadCIFAtoms(Filename,Verbose);
  if Verbose > 0
    fprintf('zReadandAnalyze: Read .cifatoms file %s\n',Filename);
  end
end

if ~isempty(ATOMNUMBER)

  if Verbose > 1
    ATOMNUMBER(1:10)
    ATOMNAME(1:10)
    VERSION(1:10)
    UNITNAME(1:10)
    CHAIN(1:10)
    NTNUMBER(1:10)
  end

  P = [P BETA];  % include beta factors with positions to make it easier to assign them to atoms

  if Readable == 0
    fprintf('zReadandAnalyze was unable to read the PDB file %s\n',Filename);
    fprintf('Please move it to FR3D/PDBFiles/Trouble Reading.\n');
  end

  % Move models in NMR file apart ---------------------------------------------
  % Otherwise they land right on top of one another and confuse basepairs
  % But dont do this with files like 3NJ6, which are x-ray and have 2 models

  if (max(ModelNum) > 2),
    P(:,1) = P(:,1) + 1000*(ModelNum - 1);  % move x by 1000 Angstroms
  end

  % Read standard bases from Jiri Sponer QM calculations ----------------------

  zStandardBases                % read in QM locations of atoms in 4 bases

  Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  base_to_atoms = containers.Map();
  base_to_atoms('A') = A_Atoms;
  base_to_atoms('C') = C_Atoms;
  base_to_atoms('G') = G_Atoms;
  base_to_atoms('U') = U_Atoms;

  sugar_atom_names = {'C1''','C2''','C3''','C4''','C5''','O2''','O3''','O4''','O5''','P','OP1','OP2','OP3'};

  % prepare a place to store modified to standard mappings, but don't load until needed
  modified_base_to_parent = containers.Map();

  % Define standard amino acid names
  AminoAcids = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','PYL','SEC','SER','THR','TRP','TYR','VAL'};  % 22 values

  % no need to print information about these or load modified nucleotides when these are encountered
  common_non_nucleotide_identifiers = {'MG','ZN','NA','HOH','UNK','K'};

  % Extract locations of base and sugar atoms ---------------------------------

  NT  = [];                                    % nucleotide data structure
  n   = 1;                                     % current NT index
  AA  = [];                                    % amino acid data structure
  aa  = 1;                                     % current amino acid index
  Het = [];                                    % HETATM data structure
  hh  = 1;                                     % current HETATM number
  i   = 1;                                     % current atom/row number
  ntDataSize  = 1000;                          % manage space allocation
  aaDataSize  = 1000;                          % manage space allocation
  hetDataSize = 1000;                          % manage space allocation

  while i <= length(NTNUMBER)                  % go through all atoms

    % accumulate HETATM entries separately, skip over them
    % that is needed for MG and other atoms
    % while i <= length(NTNUMBER) && strcmp(ATOM_TYPE{i},'HETATM')
    %   Het(hh).AtomNumber = ATOMNUMBER(i);
    %   Het(hh).Atom       = ATOMNAME{i};
    %   Het(hh).Unit       = UNITNAME{i};
    %   Het(hh).Chain      = CHAIN{i};
    %   Het(hh).Number     = NTNUMBER{i};
    %   Het(hh).Loc        = P(i,1:3);
    %   Het(hh).Beta       = P(i,4);
    %   Het(hh).Center     = P(i,1:3);
    %   Het(hh).ModelNum   = ModelNum(i);   % Anton 9/14/2011
    %   i  = i + 1;
    %   hh = hh + 1;

    %   if hh >= hetDataSize
    %     hetDataSize = hetDataSize * 2;
    %     Het(hetDataSize).Unit = '';
    %   end
    % end

    if i <= length(NTNUMBER)

      UnitType = 0;                             % RNA NT = 1, Amino acid = 2, etc.
      j = [];                                   % initialize rows of next nt
      ntnum = NTNUMBER{i};                      % nucleotide number from pdb file
      chnum = CHAIN{i};                         % Anton 12/12/10
      % to deal with consecutive nucleotides that differ only by chain:
      while (i <= length(NTNUMBER)) && (strcmp(NTNUMBER{i},ntnum)) && (strcmp(CHAIN{i},chnum))  %Anton 12/12/10
                                                % pull out rows for this unit
                                                % assumes they are continguous
        j = [j; i];                             % add rows for this nucleotide
        i = i + 1;                              % go to the next row
      end

      Sugar = Inf * ones(12,4);
      Loc   = Inf * ones(11,4);                 % make it possible to recognize missing atom locations

      UnitName = UNITNAME{j(1)};                % one, two, or three letter code
      if exist('UNITID') && length(UNITID{j(1)}) > 0
        CurrentID = UNITID{j(1)};
      else
        % make unit id; ignores insertion codes, alternate ids, symmetries
        CurrentID = sprintf('%s|%d|%s|%s|%s',PDBID,ModelNum(j(1)),CHAIN{j(1)},UNITNAME{j(1)},NTNUMBER{j(1)});
      end

      RNAbases = {'A','C','G','U','C+','+U'};   % standard RNA bases

      if ismember(UnitName,RNAbases)            % looks like an RNA nucleotide, or Y or I or something
        NT(n).Base     = UNITNAME{j(1)};        % letter for this nucleotide
        NT(n).Unit     = UNITNAME{j(1)};        % letter for this nucleotide
        NT(n).Chain    = CHAIN{j(1)};           % what chain nucleotide is in
        NT(n).Number   = NTNUMBER{j(1)};        % nucleotide number for this molecule
        NT(n).ModelNum = ModelNum(j(1));        % model number, especially for NMR
        NT(n).ID       = CurrentID;
        NT(n).AltLoc   = VERSION{j(1)};
      end

      NTBackbone = 0;                         % distinguish AAs from modified NTs

      if strcmp(UnitName,'A')
        for k = min(j):max(j)
          if strcmp(VERSION{k},'A') || strcmp(VERSION{k},'.') || strcmp(VERSION{k},'?') || isempty(VERSION{k})
            switch ATOMNAME{k}
              case 'N9',      Loc( 1,:) = P(k,:);
              case 'C4',      Loc( 2,:) = P(k,:);
              case 'N3',      Loc( 3,:) = P(k,:);
              case 'N1',      Loc( 4,:) = P(k,:);
              case 'C6',      Loc( 5,:) = P(k,:);
              case 'N6',      Loc( 6,:) = P(k,:);
              case 'C8',      Loc( 7,:) = P(k,:);
              case 'C5',      Loc( 8,:) = P(k,:);
              case 'C2',      Loc( 9,:) = P(k,:);
              case 'N7',      Loc(10,:) = P(k,:);
              case {'C1*','C1'''},     Sugar( 1,:) = P(k,:);
              case {'C2*','C2'''},     Sugar( 2,:) = P(k,:);
              case {'O2*','O2'''},     Sugar( 3,:) = P(k,:);
              case {'C3*','C3'''},     Sugar( 4,:) = P(k,:);
              case {'O3*','O3'''},     Sugar( 5,:) = P(k,:);
              case {'C4*','C4'''},     Sugar( 6,:) = P(k,:);
              case {'O4*','O4'''},     Sugar( 7,:) = P(k,:);
              case {'C5*','C5'''},     Sugar( 8,:) = P(k,:);
              case {'O5*','O5'''},     Sugar( 9,:) = P(k,:);
              case 'P',                Sugar(10,:) = P(k,:);
              case {'O1P','OP1'},      Sugar(11,:) = P(k,:);
              case {'O2P','OP2'},      Sugar(12,:) = P(k,:);
            end
          end
        end
        Loc(11,:)    = zeros(1,4);
        NT(n).Loc    = Loc(:,1:3);
        NT(n).Sugar  = Sugar(:,1:3);
        NT(n).Beta   = [Sugar(:,4); Loc(1:10,4)];   % sugar first, then base
        NT(n).Center = mean(Loc(1:10,1:3));         % mean of heavy base atoms
        NT(n).Code   = 1;                           % A is 1, C is 2, etc.
        UnitType = 1;                               % RNA nucleotide
        n = n + 1;
      elseif strcmp(UnitName,'C') || strcmp(UnitName,'C+'),
        for k = min(j):max(j),
          if strcmp(VERSION{k},'A') || strcmp(VERSION{k},'.') || strcmp(VERSION{k},'?') || isempty(VERSION{k}),
            switch ATOMNAME{k},
              case 'N1',      Loc( 1,:) = P(k,:);
              case 'C2',      Loc( 2,:) = P(k,:);
              case 'O2',      Loc( 3,:) = P(k,:);
              case 'N3',      Loc( 4,:) = P(k,:);
              case 'C4',      Loc( 5,:) = P(k,:);
              case 'N4',      Loc( 6,:) = P(k,:);
              case 'C6',      Loc( 7,:) = P(k,:);
              case 'C5',      Loc( 8,:) = P(k,:);
              case {'C1*','C1'''},     Sugar( 1,:) = P(k,:);
              case {'C2*','C2'''},     Sugar( 2,:) = P(k,:);
              case {'O2*','O2'''},     Sugar( 3,:) = P(k,:);
              case {'C3*','C3'''},     Sugar( 4,:) = P(k,:);
              case {'O3*','O3'''},     Sugar( 5,:) = P(k,:);
              case {'C4*','C4'''},     Sugar( 6,:) = P(k,:);
              case {'O4*','O4'''},     Sugar( 7,:) = P(k,:);
              case {'C5*','C5'''},     Sugar( 8,:) = P(k,:);
              case {'O5*','O5'''},     Sugar( 9,:) = P(k,:);
              case 'P',       Sugar(10,:) = P(k,:);
              case {'O1P','OP1'},     Sugar(11,:) = P(k,:);
              case {'O2P','OP2'},     Sugar(12,:) = P(k,:);
            end
          end
        end
        Loc( 9,:) = zeros(1,4);
        Loc(10,:) = zeros(1,4);
        Loc(11,:) = zeros(1,4);
        NT(n).Loc    = Loc(:,1:3);
        NT(n).Sugar  = Sugar(:,1:3);
        NT(n).Beta   = [Sugar(:,4); Loc(1:8,4)];   % sugar first, then base
        NT(n).Center = mean(Loc(1:8,1:3));
        NT(n).Code   = 2;                          % A is 1, C is 2, etc.
        UnitType = 1;                              % RNA nucleotide
        n = n + 1;
      elseif strcmp(UnitName,'G'),
        for k = min(j):max(j),
          if strcmp(VERSION{k},'A') || strcmp(VERSION{k},'.') || strcmp(VERSION{k},'?') || isempty(VERSION{k}),
            switch ATOMNAME{k},
              case 'N9',      Loc( 1,:) = P(k,:);
              case 'C4',      Loc( 2,:) = P(k,:);
              case 'N3',      Loc( 3,:) = P(k,:);
              case 'N1',      Loc( 4,:) = P(k,:);
              case 'C6',      Loc( 5,:) = P(k,:);
              case 'O6',      Loc( 6,:) = P(k,:);
              case 'C8',      Loc( 7,:) = P(k,:);
              case 'C5',      Loc( 8,:) = P(k,:);
              case 'C2',      Loc( 9,:) = P(k,:);
              case 'N7',      Loc(10,:) = P(k,:);
              case 'N2',      Loc(11,:) = P(k,:);
              case {'C1*','C1'''},     Sugar( 1,:) = P(k,:);
              case {'C2*','C2'''},     Sugar( 2,:) = P(k,:);
              case {'O2*','O2'''},     Sugar( 3,:) = P(k,:);
              case {'C3*','C3'''},     Sugar( 4,:) = P(k,:);
              case {'O3*','O3'''},     Sugar( 5,:) = P(k,:);
              case {'C4*','C4'''},     Sugar( 6,:) = P(k,:);
              case {'O4*','O4'''},     Sugar( 7,:) = P(k,:);
              case {'C5*','C5'''},     Sugar( 8,:) = P(k,:);
              case {'O5*','O5'''},     Sugar( 9,:) = P(k,:);
              case 'P',       Sugar(10,:) = P(k,:);
              case {'O1P','OP1'},     Sugar(11,:) = P(k,:);
              case {'O2P','OP2'},     Sugar(12,:) = P(k,:);
            end
          end
        end
        NT(n).Loc    = Loc(:,1:3);
        NT(n).Sugar  = Sugar(:,1:3);
        NT(n).Beta   = [Sugar(:,4); Loc(:,4)];   % sugar first, then base
        NT(n).Center = mean(Loc(1:11,1:3));
        NT(n).Code   = 3;                          % A is 1, C is 2, etc.
        UnitType = 1;                              % RNA nucleotide
        n = n + 1;
      elseif strcmp(UnitName,'U') || strcmp(UnitName,'+U'),
        for k = min(j):max(j),
          if strcmp(VERSION{k},'A') || strcmp(VERSION{k},'.') || strcmp(VERSION{k},'?') || isempty(VERSION{k}),
            switch ATOMNAME{k},
              case 'N1',      Loc( 1,:) = P(k,:);
              case 'C2',      Loc( 2,:) = P(k,:);
              case 'O2',      Loc( 3,:) = P(k,:);
              case 'N3',      Loc( 4,:) = P(k,:);
              case 'C4',      Loc( 5,:) = P(k,:);
              case 'O4',      Loc( 6,:) = P(k,:);
              case 'C6',      Loc( 7,:) = P(k,:);
              case 'C5',      Loc( 8,:) = P(k,:);
              case {'C1*','C1'''},     Sugar( 1,:) = P(k,:);
              case {'C2*','C2'''},     Sugar( 2,:) = P(k,:);
              case {'O2*','O2'''},     Sugar( 3,:) = P(k,:);
              case {'C3*','C3'''},     Sugar( 4,:) = P(k,:);
              case {'O3*','O3'''},     Sugar( 5,:) = P(k,:);
              case {'C4*','C4'''},     Sugar( 6,:) = P(k,:);
              case {'O4*','O4'''},     Sugar( 7,:) = P(k,:);
              case {'C5*','C5'''},     Sugar( 8,:) = P(k,:);
              case {'O5*','O5'''},     Sugar( 9,:) = P(k,:);
              case 'P',                Sugar(10,:) = P(k,:);
              case {'O1P','OP1'},      Sugar(11,:) = P(k,:);
              case {'O2P','OP2'},      Sugar(12,:) = P(k,:);
            end
          end
        end
        Loc( 9,:) = zeros(1,4);
        Loc(10,:) = zeros(1,4);
        Loc(11,:) = zeros(1,4);
        NT(n).Loc    = Loc(:,1:3);
        NT(n).Sugar  = Sugar(:,1:3);
        NT(n).Beta   = [Sugar(:,4); Loc(1:8,4)];   % sugar first, then base
        NT(n).Center = mean(Loc(1:8,1:3));
        NT(n).Code   = 4;                          % A is 1, C is 2, etc.
        UnitType = 1;                              % RNA nucleotide
        n = n + 1;
      elseif ismember(UnitName,AminoAcids)        % one of 22 standard amino acids
        AA(aa).AtomNumber = ATOMNUMBER(j);
        AA(aa).Atom       = ATOMNAME(j);        % store atom names
        AA(aa).Unit       = UNITNAME{j(1)};
        AA(aa).Chain      = CHAIN{j(1)};        % what chain the a.a. is in
        AA(aa).Number     = NTNUMBER{j(1)};     % number for this a.a.
        AA(aa).Loc        = P(j,1:3);
        AA(aa).Beta       = P(j,4);
        k = min(4,length(j));
        AA(aa).Center     = mean(P(j(1:k),1:3),1); % backbone center
        AA(aa).ModelNum   = ModelNum(j(1));     % Anton 9/14/2011
        AA(aa).ID         = CurrentID;

        aa = aa + 1;
        UnitType = 2;               % amino acid

        for m = 1:length(AA(aa-1).Atom)
          if ~isempty(strfind(AA(aa-1).Atom{m},''''))  % look for a prime on atom name
            NTBackbone = 1;
          end
        end

      elseif ismember(UnitName,common_non_nucleotide_identifiers)  % common ions and water
        % common ions and water, no need to load all modified nucleotides

      else
        % unrecognized unit, may be a modified nucleotide, make sure to load the definitions
        if length(modified_base_to_parent.keys()) == 0
          fprintf('zReadandAnalyze: Loading modified nucleotides\n')
          [modified_base_to_parent,modified_base_atom_list,modified_atom_to_parent,parent_atom_to_modified,modified_base_to_hydrogens,modified_base_to_hydrogen_coordinates] = zDefineModifiedNucleotides();  % load in the modified nucleotide definitions
        end

        if isKey(modified_base_to_parent, UnitName)
          fprintf('zReadandAnalyze: Storing coordinates of modified residue %s unit %s\n',CurrentID,UNITNAME{j(1)});

          standard_base = modified_base_to_parent(UnitName); % standard base for this modified NT
          modified_atom_to_parent_atom = modified_atom_to_parent(UnitName);

          % store modified nucleotide with the name of the standard nucleotide
          NT(n).Base     = standard_base;         % standard base for this modified NT
          NT(n).Unit     = UNITNAME{j(1)};        % actual name of this nucleotide
          NT(n).Chain    = CHAIN{j(1)};           % what chain nucleotide is in
          NT(n).Number   = NTNUMBER{j(1)};        % nucleotide number for this molecule
          NT(n).ModelNum = ModelNum(j(1));        % model number, especially for NMR
          NT(n).ID       = CurrentID;             % keep the original unit id, with actual name
          NT(n).AltLoc   = VERSION{j(1)};

          LocDict   = containers.Map;   % use a dictionary to store locations by atom name
          SugarDict = containers.Map;   % use a dictionary to store locations by atom name
          BetaDict  = containers.Map;   % use a dictionary to store locations by atom name

          % loop over atoms of the modified residue in the .cif file
          for k = min(j):max(j)
            if isempty(VERSION{k}) || strcmp(VERSION{k},'A') || strcmp(VERSION{k},'.') || strcmp(VERSION{k},'?')
              % do not work with alternate locations like B

              modified_atom = ATOMNAME{k};

              if isKey(modified_atom_to_parent_atom,modified_atom)
                % if there is a corresponding parent atom, store that, otherwise don't store
                parent_atom = modified_atom_to_parent_atom(modified_atom);
                if ismember(parent_atom,sugar_atom_names)
                  % does not include hydrogen atoms
                  SugarDict(parent_atom) = P(k,1:3);  % store coordinates by parent atom name
                  BetaDict(modified_atom) = P(k,4);   % store beta value by atom name
                elseif ismember(parent_atom,base_to_atoms(standard_base))
                  % includes hydrogen atoms
                  LocDict(parent_atom) = P(k,1:3);    % store coordinates by parent atom name
                  BetaDict(modified_atom) = P(k,4);   % store beta value by atom name
                end

                % try to save sugar atoms in the same way as standard bases
                switch parent_atom
                  case {'C1*','C1'''},     Sugar( 1,:) = P(k,:);
                  case {'C2*','C2'''},     Sugar( 2,:) = P(k,:);
                  case {'O2*','O2'''},     Sugar( 3,:) = P(k,:);
                  case {'C3*','C3'''},     Sugar( 4,:) = P(k,:);
                  case {'O3*','O3'''},     Sugar( 5,:) = P(k,:);
                  case {'C4*','C4'''},     Sugar( 6,:) = P(k,:);
                  case {'O4*','O4'''},     Sugar( 7,:) = P(k,:);
                  case {'C5*','C5'''},     Sugar( 8,:) = P(k,:);
                  case {'O5*','O5'''},     Sugar( 9,:) = P(k,:);
                  case 'P',                Sugar(10,:) = P(k,:);
                  case {'O1P','OP1'},      Sugar(11,:) = P(k,:);
                  case {'O2P','OP2'},      Sugar(12,:) = P(k,:);
                end
              end
            end
          end

          NT(n).Loc        = LocDict;                  % dictionary keyed by atom
          NT(n).Sugar      = Sugar(:,1:3);             % matrix of sugar atom coordinates; maybe some missing
          NT(n).SugarDict  = SugarDict;                % dictionary keyed by atom
          NT(n).Beta       = BetaDict;                 % dictionary keyed by atom
          NT(n).Center     = [0 0 0];                  % will be revised later
          NT(n).Code       = 9;                        % A is 1, C is 2, ... modified nucleotides are 9
          UnitType = 4;                                % modified RNA nucleotide
          n = n + 1;
        end
      end                                         % end A, C, G, U, amino, modified cases

      % allocate more space for NT variable when it is full
      if n >= ntDataSize
        ntDataSize = ntDataSize * 2;
        NT(ntDataSize).Unit = '';
      end

      % allocate more space for AA variable when it is full
      if aa >= aaDataSize
        aaDataSize = aaDataSize * 2;
        AA(aaDataSize).Unit = '';
      end

      if UnitType == 0 || NTBackbone == 1      % unrecognized unit
        if ismember(UNITNAME{j(1)},common_non_nucleotide_identifiers)
          % do not print anything, these are too common
        else
          fprintf('zReadandAnalyze: Unrecognized unit:  %s will be added as HETATM ========= \n',CurrentID);
        end
        if NTBackbone == 1
          aa = aa - 1;                            % not an amino acid
        end

        for hhh = 1:length(j)
          Het(hh+hhh-1).AtomNumber = ATOMNUMBER(j(hhh));
          Het(hh+hhh-1).Atom       = ATOMNAME{j(hhh)};
          Het(hh+hhh-1).Unit       = UNITNAME{j(hhh)};
          Het(hh+hhh-1).Chain      = CHAIN{j(hhh)};
          Het(hh+hhh-1).Number     = NTNUMBER{j(hhh)};
          Het(hh+hhh-1).Loc        = P(j(hhh),1:3);
          Het(hh+hhh-1).Beta       = P(j(hhh),4);
          Het(hh+hhh-1).Center     = P(j(hhh),1:3);
        end
        hh = hh + length(j);

      elseif (UnitType == 1)                     % standard RNA base

        if (max(max(Loc)) == Inf)                % base atoms missing
          n = n - 1;                              % pretend this nt was never added
          if Verbose > 0
            NumGood = length(find((Loc(1,:) < Inf) .* (abs(Loc(1,:)) > 0)));
            fprintf('Base %s has %d atoms, so it will be skipped\n',CurrentID,NumGood);
          end
        elseif (max(max(Sugar)) == Inf)          % sugar atom missing
          NumGood = 12;
          for k = 1:12
            if Sugar(k,1) == Inf
              Sugar(k,:) = Loc(1,:);              % use glycosidic atom instead
              NumGood    = NumGood - 1;
            end
          end
          NT(n-1).Sugar = Sugar(:,1:3);           % n has already been incremented, so look back to fix
          if Verbose > 0
            fprintf('Base %s backbone has %d atoms, substituting glycosidic atom\n',CurrentID,NumGood);
          end
        end
      end
    end
  end                                           % end while i < ... loop

  NumNT = n - 1;                                % total number of nucleotides
  NumAA = aa - 1;
  NumHet = hh - 1;

  % Compute best shift and rotation matrices ---------------------------------

  warning off MATLAB:divideByZero

  for n=1:NumNT                                      % analyze all nucleotides
    C  = NT(n).Code;                                 % A is 1, C is 2, ...
    % fprintf('Nucleotide %15s has code %d and unit %s\n', NT(n).ID, C, NT(n).Unit);
    if C <= 4                                        % standard RNA base
      L  = Lim(1,C);                                 % number of atoms in base
      X  = StandardLoc(1:L,:,C);                     % ideal base atom locations
      Y  = NT(n).Loc(1:L,:);                         % observed atom locations

      [r, sc, sh] = zBestTransformation(X,Y);        % find best rotation, shift

      L2 = Lim(2,C);                                 % num of atoms and hyd in base
      X2 = StandardLoc(1:L2,:,C);                    % ideal base & hyd locations
      F  = (sh*ones(1,L2) + r*X2')';                 % best fit without scaling
      e  = sqrt(sum(sum((Y - F(1:L,:)).^2)))/L;      % error measure;
                                                     % should be between 0 and 10
      if (e > 0.1) && Verbose > 0,
        fprintf('Nucleotide %c%s has average fitting error %6.4f Angstroms\n', NT(n).Base, NT(n).Number, e);
      end

      NT(n).Rot         = r;                         % save the rotation matrix
      NT(n).Fit(1:L2,:) = F;                         % fitted locations of base, H

    elseif C == 9                                    % known modified nucleotide

      modified_nucleotide = NT(n).Unit;
      % fprintf('Working on modified nucleotide %s\n', modified_nucleotide);
      parent_atom_to_modified_atom = parent_atom_to_modified(modified_nucleotide);

      C = find(NT(n).Base == 'ACGU');                % code number of parent base
      parent_base_atoms = base_to_atoms(NT(n).Base); % cell array of atom names

      X = [];                                        % parent atom locations in the plane
      Y = [];                                        % modified nt locations in space

      for a = 1:Lim(1,C)                             % loop over heavy atoms in parent base
        parent_atom = parent_base_atoms{a};          % name of parent atom
        if isKey(NT(n).Loc,parent_atom)
          X = [X; StandardLoc(a,:,C)];               % coordinates of parent atom in the plane
          Y = [Y; NT(n).Loc(parent_atom)];           % location in space of modified base atom
        end
      end

      [r, sc, sh] = zBestTransformation(X,Y);        % find best rotation, shift

      L2 = Lim(2,C);                                 % num of atoms and hyd in base
      X2 = StandardLoc(1:L2,:,C);                    % ideal base & hyd locations
      F  = (sh*ones(1,L2) + r*X2')';                 % map standard base atoms into 3D space,
                                                     % even if the modified base does not have them
%      e  = sqrt(sum(sum((Y - F(1:L,:)).^2)))/L;      % error measure; too hard for mod nt

      NT(n).Rot         = r;                         % save the rotation matrix
      NT(n).Fit(1:L2,:) = F;                         % fitted parent base atoms
      NT(n).Center      = mean(F(1:Lim(1,C),:));     % use center of fitted parent heavy atoms

      if Verbose > 1
        figure(1)
        clf
        VP.LabelBases = 10;
        VP.LabelAtoms = 10;
        VP.Sugar = 1;
        zPlotOneNT(NT(n),VP)
        pause
      end
    end
    NT(n).Syn         = 0;                           % determine anti/syn later
  end

  % Fill in fields of File ----------------------------------------------------

  File.Filename = CoordinateFile;
  File.Filename  = Filename;
  File.NT        = NT(1:NumNT);
  File.NumNT     = NumNT;
  File.Distance  = [];
  File.HandClass = [];
  File.Comment   = [];
  File.CI        = sparse(NumNT,NumNT);
  File.Edge      = sparse(NumNT,NumNT);
  File.Modified  = 1;
  File.Header    = [];
  File.BasePhosphate = [];
  File.AA        = AA(1:NumAA);
  File.Het       = Het(1:NumHet);

  % Calculate configuration (syn or anti) -------------------------------------

  if ~isempty(File.NT),
    SynList = mSynList(File);

    j = find(SynList);

    for k=1:length(j),
      File.NT(j(k)).Syn = 1;
    end
  end

else

  fprintf('zReadandAnalyze: Could not open file %s\n', Filename);
  NumNT = 0;
  File.Filename = Filename;
  File.NT        = [];
  File.NumNT     = -1;
  File.Distance  = [];
  File.HandClass = [];
  File.Comment   = [];
  File.CI        = sparse(NumNT,NumNT);
  File.Edge      = sparse(NumNT,NumNT);
  File.Coplanar  = sparse(NumNT,NumNT);
  File.Modified  = 0;
  File.Pair      = [];
  File.ClassVersion = 0;
  File.Header    = [];
  File.Info.Resolution  = [];
  File.Info.Descriptor  = '';
  File.Info.ExpTechnique= '';
  File.Info.ReleaseDate = '';
  File.Info.Author      = '';
  File.Info.Keywords    = '';
  File.Info.Source      = '';
  File.BasePhosphate = sparse(NumNT,NumNT);
  File.AA               = [];
  File.Het              = [];
end

File = orderfields(File);
