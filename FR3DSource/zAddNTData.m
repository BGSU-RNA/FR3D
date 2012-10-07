% zAddNTData(Filenames,ReadCode,File,Verbose,PDBStart) reads RNA structure 
% data files, if
% necessary, so that all molecules listed in Filenames are present in File.
% The parameter Index is a non-redundant list of indices of File
% corresponding to names in Filenames.

% ReadCode = 0,1,2 : load .mat files
% ReadCode = 3     : load, but do not append to File (for reclassification)
% ReadCode = 4     : load _small.mat files, do not compute distances

% F = zAddNTData('AllFiles_list',3,[],1);
% F = zAddNTData('AllFiles_list',3,[],1,'back');
% F = zAddNTData('NonRedundant_2008_02_21_list',0,[],1);

function [File,Index] = zAddNTData(Filenames,ReadCode,File,Verbose,PDBStart)

if nargin < 2,
  ReadCode = 0;                           % default is to read .mat files
end

if nargin < 4,
  Verbose = 0;
end

LoadedFiles = {};
F = 0;

if nargin >= 3,
  F = length(File);
  for j = 1:length(File),
    LoadedFiles{j} = lower(File(j).Filename);
    if isempty(LoadedFiles{j}),
      LoadedFiles{j} = '';                     % use empty string
    end
  end
end

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};                % make into a cell array
end

% ----------------------------------------- Read PDB lists, if any

FullList = [];

for j=1:length(Filenames),
  FullList = [FullList; zReadPDBList(Filenames{j},1)];
end

% ----------------------------------------- Skip some files

if nargin == 5,
  if strcmp(PDBStart,'back') == 1,
    FullList = FullList(end:-1:1);
  else
    keep = [];
    for j=1:length(FullList),
      if issorted([lower(PDBStart(1:4)); lower(FullList{j}(1:4))],'rows'),
        keep = [keep j];
      end
    end
    FullList = FullList(keep);
  end
end

% ----------------------------------------- Read PDB files

if length(FullList) > 0,

for f = 1:length(FullList),                       % loop through PDB list
  if ~isempty(FullList{f}),
    i = strmatch(lower(FullList{f}), LoadedFiles, 'exact');
    if isempty(i),                                  % if PDB not loaded,
      NewF = zGetNTData(FullList{f},ReadCode,Verbose); %   load it
      if ReadCode ~= 3,
        if F == 0,
          clear File
          File(1) = NewF;
        else
          File(F+1) = NewF;
        end
      end
      clear NewF;
      F = length(File);
      k = length(LoadedFiles);
      LoadedFiles{k+1} = FullList{f};
      Index(f) = F;                           %   point to it
    else                                      % but if PDB has been loaded
      Index(f) = i(1);                        %   point to first instance
      if length(File(i(1)).NT) == 0,
        NewF = zGetNTData(File(Index(f)).Filename,ReadCode,Verbose);
        if ReadCode ~= 3,
          File(Index(f)) = NewF;
          clear NewF;
        end
      end
    end
  end
end

% -----------------------------------  allow for File(Index)
F = length(File);

for i = 1:length(Index),             
  if Index(i) == 0,
    File(F+1).Filename = 'Fictitious';
    File(F+1).NumNT = 0;               % create a fictitious file
    Index(i) = F+1;                    % point to the fictitious file
  end
end

else

fprintf('No files specified to read in %s\n', Filenames{1});

end
