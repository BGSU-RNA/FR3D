ModelNum = [];

if isfield(Header,'Expdata') & (length(Header.ModelStart) > 1)
  if ~isempty(strfind(Header.Expdata,'NMR')),
    Header.ModelStart = [Header.ModelStart length(NTNUMBER)+1];
    for m = 1:(length(Header.ModelStart)-1),
      a = Header.ModelStart(m);
      b = Header.ModelStart(m+1)-1;
      P(a:b,1) = P(a:b,1) + (m-1)*1000;     % shift by 1000 for each model
      ModelNum(a:b,1) = m;                  % record model number
    end
  end
end

% Read standard bases from Sponer QM calculations --------------------------

zStandardBases                % read in QM locations of atoms in 4 bases

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

% Extract locations of base and sugar atoms ---------------------------------

NT  = [];                                    % nucleotide data structure
n   = 1;                                     % current NT index
AA  = [];                                    % amino acid data structure
Het = [];                                    % HETATM data structure
aa  = 1;                                     % current amino acid index
i   = 1;                                     % current atom/row number
hh  = 1;                                     % current HETATM number

while i <= length(NTNUMBER),                 % go through all atoms

 while i <= length(NTNUMBER) && strcmp(ATOM_TYPE{i},'HETATM'),
   Het(hh).AtomNumber = ATOMNUMBER(i);
   Het(hh).Atom       = ATOMNAME{i};
   Het(hh).Unit       = NTLETTER{i};
   Het(hh).Chain      = CHAIN{i};
   Het(hh).Number     = NTNUMBER{i};
   Het(hh).Loc        = P(i,1:3);
   Het(hh).Beta       = P(i,4);
   Het(hh).Center     = P(i,1:3);
   i  = i + 1;
   hh = hh + 1;
 end

 if i < length(NTNUMBER),

  UnitType = 0;                             % RNA NT = 1, Amino acid = 2, etc.
  j = [];                                   % initialize rows of next nucleo.
  ntnum = NTNUMBER{i};                      % nucleotide number from pdb file

  while (i <= length(NTNUMBER)) & (strcmp(NTNUMBER{i},ntnum)), 
                                            % pull out rows for this nucleotide
                                            % assumes they are continguous
    j = [j; i];                             % add rows for this nucleotide
    i = i + 1;                              % go to the next row
  end

  Sugar = Inf * ones(12,4);
  Loc   = Inf * ones(11,4);

  NTBase = NTLETTER{j(1)};                % one or three letter code

  if isempty(ModelNum),
    NT(n).Model = 1;                      % default model number, not NMR
  else
    NT(n).Model = ModelNum(j(1));         % use current model number for NMR
  end

  if length(NTBase) == 1,                 % looks like a nucleotide
    NT(n).Base    = NTLETTER{j(1)};       % letter for this nucleotide
    NT(n).Unit    = NTLETTER{j(1)};       % letter for this nucleotide
    NT(n).Chain   = CHAIN{j(1)};          % what chain nucleotide is in
    NT(n).Number  = NTNUMBER{j(1)};       % nucleotide number for this molecule
  end

  NTBackbone = 0;                         % distinguish AAs from modified NTs

  if strcmp(NTBase,'A'),
    for k = min(j):max(j),
     if strcmp(VERSION{k},'A') || isempty(VERSION{k}),
      switch ATOMNAME{k},
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
    UnitType = 1;                              % RNA nucleotide
    n = n + 1;
  elseif strcmp(NTBase,'C') || strcmp(NTBase,'C+'),
    for k = min(j):max(j),
     if strcmp(VERSION{k},'A') || isempty(VERSION{k}),
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
  elseif strcmp(NTBase,'G'),
    for k = min(j):max(j),
     if strcmp(VERSION{k},'A') || isempty(VERSION{k}),
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
  elseif strcmp(NTBase,'U') || strcmp(NTBase,'+U'),
    for k = min(j):max(j),
     if strcmp(VERSION{k},'A') || isempty(VERSION{k}),
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
  elseif length(NTBase) == 3,                  % probably an amino acid

    AA(aa).AtomNumber = ATOMNUMBER(j);
    AA(aa).Atom       = ATOMNAME(j);        % store atom names
    AA(aa).Unit       = NTLETTER{j(1)};
    AA(aa).Chain      = CHAIN{j(1)};        % what chain the a.a. is in
    AA(aa).Number     = NTNUMBER{j(1)};     % number for this a.a.
    AA(aa).Loc        = P(j,1:3);
    AA(aa).Beta       = P(j,4);
    k = min(4,length(j));
    AA(aa).Center     = mean(P(j(1:k),1:3)); % backbone center

    aa = aa + 1;
    UnitType = 2;

    % ------------------------- Might this be a nucleotide with 3 letter name?

    for m = 1:length(AA(aa-1).Atom),
      if ~isempty(strfind(AA(aa-1).Atom{m},'''')),
        NTBackbone = 1;
      end
    end
  end                                         % end four if statements

  if (UnitType == 0) || NTBackbone == 1,      % unrecognized nucleotide
    fprintf('Unrecognized unit:  %s%4s_%s will be added as HETATM ========= \n',NTBase,NTNUMBER{j(1)},CHAIN{j(1)});
    if NTBackbone == 1,
      aa = aa - 1;                            % not an amino acid
    end

    for hhh = 1:length(j),
      Het(hh+hhh-1).AtomNumber = ATOMNUMBER(j(hhh));
      Het(hh+hhh-1).Atom       = ATOMNAME{j(hhh)};
      Het(hh+hhh-1).Unit       = NTLETTER{j(hhh)};
      Het(hh+hhh-1).Chain      = CHAIN{j(hhh)};
      Het(hh+hhh-1).Number     = NTNUMBER{j(hhh)};
      Het(hh+hhh-1).Loc        = P(j(hhh),1:3);
      Het(hh+hhh-1).Beta       = P(j(hhh),4);
      Het(hh+hhh-1).Center     = P(j(hhh),1:3);
    end
    hh = hh + length(j);

  elseif (UnitType == 1),                     % standard RNA base

    if (max(max(Loc)) == Inf),                 % base atoms missing
      n = n - 1;
      if Verbose > 0,
        NumGood = length(find((Loc(1,:) < Inf) .* (abs(Loc(1,:)) > 0)));
        fprintf('Base %s%4s_%s has %d atoms, so it will be skipped\n',NTBase,NTNUMBER{j(1)},CHAIN{j(1)},NumGood);
      end
    elseif (max(max(Sugar)) == Inf),          % sugar atom missing
      NumGood = 12;
      for k = 1:12,
        if Sugar(k,1) == Inf,
          Sugar(k,:) = Loc(1,:);              % use glycosidic atom
          NumGood    = NumGood - 1;
        end
      end
      NT(n-1).Sugar = Sugar(:,1:3);
      if Verbose > 0,
        fprintf('Base %s%4s_%s backbone has %d atoms, substituting glycosidic atom\n',NT(n-1).Base,NT(n-1).Number,NT(n-1).Chain,NumGood);
      end
    end
  end
 end
end                                           % end while i < ... loop

NumNT = n - 1;                                % total number of nucleotides
NumAA = aa - 1;

% Compute best shift and rotation matrices ---------------------------------

warning off MATLAB:divideByZero

for n=1:NumNT,                                   % analyze all nucleotides
  C  = NT(n).Code;                               % A is 1, C is 2, ...
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
  NT(n).Syn         = 0;
end

% Fill in fields of File ----------------------------------------------------

File.PDBFilename = PDBFilename;
File.Filename  = Filename;
File.NT        = NT(1:NumNT);
File.NumNT     = NumNT;
File.Distance  = [];
File.HandClass = [];
File.Comment   = [];
File.CI        = sparse(NumNT,NumNT);
File.Edge      = sparse(NumNT,NumNT);
File.Modified  = 1;
File.Header    = Header;
File.BasePhosphate = [];
File.AA        = AA;
File.Het       = Het;

% Calculate configuration (syn or anti) -------------------------------------

if length(File.NT) > 0,
  SynList = mSynList(File);

  j = find(SynList);

  for k=1:length(j),
    File.NT(j(k)).Syn = 1;
  end
end


File = orderfields(File);

