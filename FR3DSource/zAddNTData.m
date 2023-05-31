% [File,Index] = zAddNTData(Filenames,ReadCode,File,Verbose,KeepAA) reads RNA structure
% data files, if necessary, so that all molecules listed in Filenames are present in File.
% The parameter Index is a non-redundant list of indices of File
% corresponding to names in Filenames.
% When a file is read for the first time, a binary file is saved with the same filename and extension .mat.
% The default is to look for Filename in the PrecomputedData folder in the current directory.
% If not found there, the default is to look for a .cif file in the PDBFiles folder in the current directory,
% and the resulting binary .mat file is stored in the PrecomputedData folder in the current directory.
% If a full path is given as part of Filename, that path will be searched first for the .mat file, then .cif,
% and the resulting .mat file will be written to that path with the same filename.
% If Filename includes the extension .pdb or .cif, that file is used instead of the .mat file.

% ReadCode = 0 : load .mat files
% ReadCode = 1 : load .mat files and reclassify
% ReadCode = 3 : load, but do not append to File (for reclassification)
% ReadCode = 4 : read .cifatoms or .cif or .pdb file and reclassify

% F = zAddNTData('NonRedundant_2008_02_21_list',0,[],1);
% F = zAddNTData('http://rna.bgsu.edu/rna3dhub/nrlist/download/2.108/4.0A/csv',0,[],1,0)  % download representative set and use it
% F = zAddNTData('http://rna.bgsu.edu/rna3dhub/nrlist/download/current/4.0A/csv',0,[],1,0)  % download representative set and use it

% F = zAddNTData('1S72')
% F = zAddNTData('C:\Users\zirbel\Documents\FR3D\PrecomputedData\1S72')

% KeepAA is 1 one by default, but when set to 0, the amino acid field is removed from each file, reducing memory usage significantly

function [File,Index] = zAddNTData(Filenames,ReadCode,File,Verbose,KeepAA)

if nargin < 2,
    ReadCode = 0;                           % default is to read .mat files
end

if nargin < 4,
    Verbose = 0;
end

if nargin < 5,
    KeepAA = 1;
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

FullList = {};

for j=1:length(Filenames),
    FullList = [FullList zReadPDBList(Filenames{j},1)];
end

% ----------------------------------------- Read PDB files

if length(FullList) > 0,

    for f = 1:length(FullList),                       % loop through PDB list
        if ~isempty(FullList{f}),
            i = strmatch(lower(FullList{f}), LoadedFiles, 'exact');
            if isempty(i),                                  % if PDB not loaded,
                NewF = zLoadIFE(FullList{f},ReadCode,Verbose); % this will call zGetNTData

                if KeepAA == 0,
                    NewF = rmfield(NewF,'AA');
                end
                if ReadCode ~= 3,
                    if F == 0,
                        clear File
                        File(1) = NewF;
                    else
                        try
                            File(F+1) = NewF;
                        catch
                            fprintf('This file differs from the others and could not be added to the rest (file %d):\n',f)
                            File
                            NewF
                        end
                    end
                end
                clear NewF;
                F = length(File);
                Index(f) = F;                           % point to it
                k = length(LoadedFiles);
                LoadedFiles{k+1} = FullList{f};
            else                                        % but if PDB has been loaded
                Index(f) = i(1);                        % point to first instance
                if length(File(i(1)).NT) == 0,          % no nucleotides in the file for some reason
                    NewF = zLoadIFE(File(Index(f)).Filename,ReadCode,Verbose);
                    if KeepAA == 0,
                        NewF = rmfield(NewF,'AA');
                    end
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
