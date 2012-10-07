% zReadandAnalyze.m read a _mod.pdb file, then
% pulls out the locations of base atoms,
% and computes the best shift and rotation from a standard base

% It returns a variable File with fields:
%
% NT        array of nucleotide structures
% NumNT     the total number of distinct nucleotides found
% Distance  sparse matrix of center-center distances, if less than 15

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

function [File] = zReadandAnalyze(PDBFilename,Verbose)

i = strfind(PDBFilename,'.');
Filename = PDBFilename(1:(i(end)-1));

if nargin < 2,
  Verbose = 1;
end

% read PDBFilename ----------------------------------------------------

if exist([pwd filesep 'PDBFiles' filesep 'Trouble reading']) == 7,
  if exist([pwd filesep 'PDBFiles' filesep 'Trouble reading' filesep PDBFilename]) == 2,
    Stop = 1;
    if Verbose > 0,
      fprintf('zReadandAnalyze is skipping %s because it appeared in FR3D/PDBFiles/Trouble Reading\n', PDBFilename);
    end
  else
    Stop = 0;
  end
else
  Stop = 0;
end

fid = fopen(PDBFilename,'r');

if fid < 0,
  NewPDBFilename = strrep(PDBFilename,'.pdb','.ent');
  NewPDBFilename = ['pdb' NewPDBFilename];

  fid = fopen(NewPDBFilename,'r');

  if fid > 0,
    PDBFilename = NewPDBFilename;
  end
end


if (fid > 0) && (Stop == 0),

fclose(fid);

a = 1000+round(8900*rand);

TempFileName = ['##TempPDB' num2str(a) '.pdb'];

Header = zExtractAtomsPDB(PDBFilename,TempFileName,Verbose);

[ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, NTLETTER, CHAIN, NTNUMBER, P, Readable] = zReadPDBTextRead(TempFileName);

delete(TempFileName);

if Readable == 0,
  fprintf('zReadandAnalyze was unable to read the PDB file %s\n',PDBFilename);
  fprintf('Please move it to FR3D/PDBFiles/Trouble Reading.\n');
end

% Move models in NMR file apart ---------------------------------------------

if isfield(Header,'Expdata') & (length(Header.ModelStart) > 1)
  if ~isempty(strfind(Header.Expdata,'NMR')),
    Header.ModelStart = [Header.ModelStart length(NTNUMBER)+1];
    for m = 1:(length(Header.ModelStart)-1),
      a = Header.ModelStart(m);
      b = Header.ModelStart(m+1)-1;
      P(a:b,1) = P(a:b,1) + (m-1)*1000;     % shift by 1000 for each model
    end
  end
end

% Read standard bases from Sponer QM calculations --------------------------

zStandardBases                % read in QM locations of atoms in 4 bases

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

% Extract locations of base and sugar atoms ---------------------------------

NT = [];                                    % nucleotide data structure
n  = 1;                                     % current NT index
i  = 1;                                     % current atom/row number
unrec = 0;                                  % count unrecognized nucleotides

while i < length(NTNUMBER),                 % go through all atoms
  Flag = 0;                                 % not a recognized nucleotide
  j = [];                                   % initialize rows of next nucleo.
  ntnum = NTNUMBER{i};                      % nucleotide number from pdb file

  while (i <= length(NTNUMBER)) & (strcmp(NTNUMBER{i},ntnum)), 
                                            % pull out rows for this nucleotide
                                            % assumes they are continguous
    j = [j; i];                             % add rows for this nucleotide
    i = i + 1;                              % go to the next row
  end

  NT(n).Base    = NTLETTER{j(1)};         % letter for this nucleotide
  NT(n).Chain   = CHAIN{j(1)};            % what chain nucleotide is in
  NT(n).Number  = NTNUMBER{j(1)};         % nucleotide number for this molecule

  Sugar = Inf * ones(12,3);
  Loc   = Inf * ones(11,3);

  if (NT(n).Base == 'A') | (NT(n).Base == 'A'),
    for k = min(j):max(j),
     if (VERSION{k}=='A') | (VERSION{k}==' '),
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
        case 'P',       Sugar(10,:) = P(k,:);
        case {'O1P','OP1'},     Sugar(11,:) = P(k,:);
        case {'O2P','OP2'},     Sugar(12,:) = P(k,:);
      end
     end
    end
    Loc(11,:)    = zeros(1,3);
    NT(n).Loc    = Loc;
    NT(n).Sugar  = Sugar;
    NT(n).Center = mean(Loc(1:10,:));          % mean of heavy base atoms
    NT(n).Code   = 1;                          % A is 1, C is 2, etc.
  elseif (strcmp(NT(n).Base,'C')) | (strcmp(NT(n).Base,'C+')),
    for k = min(j):max(j),
     if (VERSION{k}=='A') | (VERSION{k}==' '),
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
    Loc( 9,:) = zeros(1,3);
    Loc(10,:) = zeros(1,3);
    Loc(11,:) = zeros(1,3);
    NT(n).Loc    = Loc;
    NT(n).Sugar  = Sugar;
    NT(n).Center = mean(Loc(1:8,:));
    NT(n).Code   = 2;                          % A is 1, C is 2, etc.
  elseif (NT(n).Base == 'G') | (NT(n).Base == 'G'),
    for k = min(j):max(j),
     if (VERSION{k}=='A') | (VERSION{k}==' '),
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
    NT(n).Loc    = Loc;
    NT(n).Sugar  = Sugar;
    NT(n).Center = mean(Loc(1:11,:));
    NT(n).Code   = 3;                          % A is 1, C is 2, etc.
  elseif (NT(n).Base == 'U') | (strcmp(NT(n).Base,'+U')),
    for k = min(j):max(j),
     if (VERSION{k}=='A') | (VERSION{k}==' '),
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
        case 'P',       Sugar(10,:) = P(k,:);
        case {'O1P','OP1'},     Sugar(11,:) = P(k,:);
        case {'O2P','OP2'},     Sugar(12,:) = P(k,:);
      end
     end
    end
    Loc( 9,:) = zeros(1,3);
    Loc(10,:) = zeros(1,3);
    Loc(11,:) = zeros(1,3);
    NT(n).Loc    = Loc;
    NT(n).Sugar  = Sugar;
    NT(n).Center = mean(Loc(1:8,:));
    NT(n).Code   = 4;                          % A is 1, C is 2, etc.
  else 
    if unrec < 5,
      if Verbose > 0,
        fprintf('Unrecognized: %s %s\n', NT(n).Base,NT(n).Number);
      end
      unrec = unrec + 1;
    end
    n = n - 1;                                 % not a recognized
                                               % nucleotide, do nothing
    Flag = 1;
  end                                         % end four if statements

  if (Flag < 1),
    if (max(max(Loc)) == Inf),                 % base atoms missing
      if Verbose > 0,
        NumGood = length(find((Loc(1,:) < Inf) .* (abs(Loc(1,:)) > 0)));
        fprintf('Base %s%s has %d atoms, so it will be skipped\n',NT(n).Base,NT(n).Number,NumGood);
      end
      n = n - 1;
    elseif (max(max(Sugar)) == Inf),          % sugar atom missing
      for k = 1:12,
        if Sugar(k,1) == Inf,
          Sugar(k,:) = Loc(1,:);              % use glycosidic atom
        end
      end
      NT(n).Sugar = Sugar;
    end
  end

  n = n + 1;                                  % next nucleotide index

end                                           % end while i < ... loop

NumNT = n - 1;                                % total number of nucleotides

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

% Calculate configuration (syn or anti) -------------------------------------

if length(File.NT) > 0,
  SynList = mSynList(File);

  j = find(SynList);

  for k=1:length(j),
    File.NT(j(k)).Syn = 1;
  end
end

else

  fprintf('Could not open file %s\n', PDBFilename);
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

end

File = orderfields(File);

end
