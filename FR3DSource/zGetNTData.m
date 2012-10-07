% zGetNTData(Filenames,ReadCode,Verbose) reads data files, depending on ReadCode
%
% If ReadCode = 0, it looks for Filename.mat and reads it if it exists.
% If ReadCode = 1, it looks for Filename.mat, reads it, and re-does 
%    the classification of interacting pairs
% If ReadCode = 2, it looks for Filename.mat, reads it 
% If ReadCode = 3, it looks for Filename.mat, reads it, re-does the
%    classification of pairs
% If ReadCode = 4, it reads Filename.pdb, analyzes each nucleotide, 
%    and classifies interacting pairs

function [Files] = zGetNTData(Filenames,ReadCode,Verbose)

[CL, CurrentVersion] = zClassLimits;  % read current version number

if nargin < 2,
  ReadCode = 0;
end

if nargin < 3,
  Verbose = 0;
end

% path(path,pwd);

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

for f=1:length(Filenames),
  Filename = Filenames{f};
  FILENAME = upper(Filename);
  filename = lower(Filename);

  ClassifyCode = 0;
  SaveCode     = 0;

  if (ReadCode > 0),
    ReadFull = 1;
  else
    ReadFull = 0;
  end

  clear File

  if ReadCode == 4,                     % re-read the PDB file
    File = zReadandAnalyze(Filename,Verbose);   % might not work on a Mac
    ClassifyCode = 1;
    ReadFull = 0;                       % no need to read PDB file again
  else                                  % try to load a precomputed version
    if (exist(strcat(Filename,'.mat'),'file') > 0),
      load(strcat(Filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded %s\n', [Filename '.mat']);
      end
    elseif (exist(strcat(FILENAME,'.MAT'),'file') > 0),  % helps on a Mac
      load(strcat(FILENAME,'.MAT'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded  %s\n', [FILENAME '.MAT']);
      end
    elseif (exist(strcat(filename,'.mat'),'file') > 0),  % helps on a Mac
      load(strcat(filename,'.mat'),'File','-mat');
      if Verbose > 0,
        fprintf('Loaded   %s\n', [filename '.MAT']);
      end
    else
      ReadFull = 1;
    end
  end

  if (ReadFull == 1),                     % try to load the full version
    File = zReadandAnalyze(Filename,Verbose);
    ClassifyCode = 1;
  end

  if ~exist('File'),
    File = [];
  end

  if ~isfield(File,'Info'),
    File = zGetPDBInfo(File);          % get resolution and other info
    SaveCode = 1;
  else
    if isempty(File.Info.Descriptor),
      File = zGetPDBInfo(File);          % look for file information
      if ~isempty(File.Info.Descriptor),
        SaveCode = 1;
      end
    end
  end

  if ~isfield(File,'ClassVersion'),
    File.ClassVersion = 0;
  end

  if ~isfield(File,'BasePhosphate'),
    File.BasePhosphate = 0;
  end

  if isfield(File,'Inter'),
    File = rmfield(File,'Inter');
  end

  if isfield(File,'SizeCode'),
    File = rmfield(File,'SizeCode');
  end

  if isfield(File,'Pair'),
    File = rmfield(File,'Pair');
  end

  if isfield(File,'Header'),
    File = rmfield(File,'Header');
  end

  Overlap = 0;

  File.Distance = [];                        % only compute this if needed

  if length(File.NT) > 1,                    % if it has nucleotides,
    if length(File.NT(1).Sugar(:,1)) < 13,
      File = zStoreO3(File);
    end

    if (ReadCode == 1) | (ReadCode == 3) | (ReadCode == 4) | ... 
      (ClassifyCode == 1) | (File.ClassVersion < CurrentVersion),

      File.Edge = sparse(File.NumNT,File.NumNT);

      c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
      File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
      d = sort(nonzeros(File.Distance));

      if ~isempty(d),

       if d(min(10,length(d))) < 1,
         fprintf('%s has overlapping nucleotides and should be avoided as such\n',File.Filename);
         Overlap = 1;
       else
         t = cputime;
         File = zClassifyPairs(File,Verbose);

         if Verbose > 0,
           fprintf(' Base-phosphate interactions ...');
         end

         File = zPhosphateInteractions(File);
         File.ClassVersion = CurrentVersion;
         ClassifyCode = 1;

         if Verbose > 0,
           fprintf(' Interaction range ... ');
         end

         File = zInteractionRange(File,Verbose);

         if Verbose > 0,
           fprintf('\nClassification took %4.2f minutes\n', (cputime-t)/60);
         end
       end
      end
    end

    if ~isfield(File.NT(1),'Syn'),
      SynList = mSynList(File);
      for k=1:length(File.NT),
        File.NT(k).Syn = SynList(k);
      end
      ClassifyCode = 1;
    end
  else
    File.ClassVersion = CurrentVersion;
  end

  if ~isfield(File,'Range'),
    File = zInteractionRange(File,Verbose);
    ClassifyCode = 1;
  end

  File = orderfields(File);

  Saved = 0;

  if (ReadCode > 0) || (ClassifyCode > 0) || (SaveCode > 0),
    zSaveNTData(File,Verbose);
    Saved = 1;
  end

  Files(f) = File;

end
