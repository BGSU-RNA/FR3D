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

if exist('OCTAVE_VERSION', 'builtin') == 5,
  warning('off');
end

% path(path,pwd);

if strcmp(class(Filenames),'char'),
  Filenames = {Filenames};
end

for f=1:length(Filenames),
  Filename = Filenames{f};

  [pathstr, name, extension] = fileparts(Filename);

  PDBID = name;

  MatFile = [PDBID '.mat'];

  ClassifyCode = 0;
  SaveCode     = 0;
  ReadText     = 0;

  if (ReadCode > 3),
    ReadFull = 1;
  else
    ReadFull = 0;
  end

  clear File

  if ~isempty(strfind(upper(extension),'.PDB')) || ~isempty(strfind(upper(extension),'.CIF')) || ~isempty(strfind(upper(extension),'.CIFATOMS')), 
    File = zReadandAnalyze(Filename);
    ReadCode = 4;
  else

    if ReadCode == 4,
      File = zReadandAnalyze(Filename,Verbose);               % read text file; trust that its filename is correct
      File.Filename = Filename;
      ClassifyCode = 1;                                       % need to classify interactions
      ReadText = 1;                                           % read the text file of coordinates
    end

    if ReadCode < 4,
      if length(pathstr) == 0,                  % try to load a .mat file; faster
        try
          load(MatFile,'File','-mat');                          % read from Matlab path
        catch
          try
            load([lower(PDBID) '.mat'],'File','-mat');
          catch
            try
              load(PDBID,'File','-mat');
            catch
              ReadCode = 5;                                     % failed to load a mat file, now look for a text file to read
            end
          end
        end
      else
        try
          load(Filename,'File','-mat');
        catch
          fprintf('zGetNTData: Unable to load %s\n',Filename);
        end
      end
    end

    if ReadCode < 4,
      if Verbose > 0,
        fprintf('Loaded  %s\n', [Filename '.mat']);
        zFlushOutput;
      end
    end

    if ReadCode == 5,
      try
        PDBFilename = [PDBID '.cifatoms'];                      % try to read the processed file with unit ids added and symmetry operations done
        File = zReadandAnalyze(PDBFilename,Verbose);
      catch
        File = zReadandAnalyze(Filename,Verbose);
      end
      File.Filename = PDBID;
      ClassifyCode = 1;                                       % need to classify interactions
      ReadText = 1;                                           % read the text file of coordinates
    end
  end

  if ~exist('File'),
    File = [];
  end

  File.PDBID = PDBID;

  if ~isfield(File,'ClassVersion'),
    File.ClassVersion = 0;
  end

  blank = sparse([],[],[],length(File.NT),length(File.NT));

  if ~isfield(File,'BasePhosphate'),
    ClassifyCode = 1;
    File.BasePhosphate = blank;
  end

  if ~isfield(File,'BaseRibose'),
    ClassifyCode = 1;
  end

  if ~isfield(File,'Coplanar'),
    ClassifyCode = 1;
    File.Coplanar = blank;
  end

  if ~isfield(File,'Covalent'),
    if length(File.NT) > 1,
      File = zBackboneContinuity(File);
    else
      File.Covalent = blank;
    end
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

  if length(File.NT) > 0,
    if ~isfield(File.NT,'Beta'),
      File.NT(1).Beta = [];            % Beta factors might not have been read
    end
  end

  if ~isfield(File,'AA'),
    File.AA = [];
  end

  if ~isfield(File,'Het'),
    File.Het = [];
  end

  Overlap = 0;

  File.Distance = [];                        % only compute this if needed

  if length(File.NT) > 1,                    % if it has nucleotides,
    if length(File.NT(1).Sugar(:,1)) < 13,
      File = zStoreO3(File);
    end

    if (ReadCode == 1) || (ReadCode == 3) || (ReadCode == 4) || ...
      (ClassifyCode == 1) || (File.ClassVersion < CurrentVersion),

      blank2 = sparse([],[],[],File.NumNT,File.NumNT);

      File.Edge = blank2;
      File.Coplanar = blank2;

      c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
      File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
      File.Distance(1,1) = 10;                 % make sure one is non-zero
      d = sort(nonzeros(File.Distance));

      if isempty(d),
       % more than one nucleotide, but too far apart to pair
       File.Range = blank2;
       File.BasePhosphate = blank2;
       File.BaseRibose = blank2;
      else
       if d(min(10,length(d))) < 1,
         fprintf('%s has overlapping nucleotides and should be used with caution\n',File.Filename);

         Overlap = 1;
         File.Range = blank2;
         File.BasePhosphate = blank2;
         File.BaseRibose = blank2;

       else

         t = cputime;
         File = zClassifyPairs(File,Verbose);

         File = zEdgeMakesMultiplePairs(File,0);

         if Verbose > 0,
           fprintf(' Base-backbone interactions ...');
         end

         File = zPhosphateInteractions(File);

         if Verbose > 0,
           fprintf(' Base-phosphate interactions ...');
         end

         File = zBaseRiboseInteractions(File);
         ClassifyCode = 1;

         if Verbose > 0,
           fprintf(' Crossing number ... ');
         end

         File = zInteractionRange(File,Verbose);

         if Verbose > 0,
           fprintf('\nClassification took %4.2f minutes\n', (cputime-t)/60);
           zFlushOutput;
         end
        end
      end
    end
  end

  File.ClassVersion = CurrentVersion;

  if length(File.NT) > 0,
   if ~isfield(File.NT(1),'glycosidicbondorientation'),
    File = zSynAnti(File);
    SaveCode = 1;
   end
  end

  if ~isfield(File,'PDBFilename'),
    File.PDBFilename = [File.Filename '.pdb'];   % will self-correct the
                                                 % next time PDB is read
  end

  if ~isfield(File,'BaseRibose'),
    File = zBaseRiboseInteractions(File);
    SaveCode = 1;
  end

  if ~isfield(File,'Backbone'),
    File = zBackboneConformation(File,Verbose);
    if File.NumNT > 1,
      File.Backbone(1,1) = 1;                    % register that we checked
    end
    SaveCode = 1;
  end

  if File.NumNT > 1 && sum(sum(File.Backbone)) == 0,
    File = zBackboneConformation(File,Verbose);
    SaveCode = 1;
  end

  if ~isfield(File,'Range') || ~isfield(File,'Crossing') || ClassifyCode > 0,
    File = zInteractionRange(File,Verbose);
    ClassifyCode = 1;
  end

  if ~isfield(File,'Flank') || ClassifyCode > 0,
    File = zBorderSS(File);                      % identify single-stranded regions
  end

  if ~isfield(File,'Redundant') || ~isfield(File,'LongestChain') || ClassifyCode > 0,
    File = zMarkRedundantNucleotides(File,Verbose);   % mark redundant and determine longest chain
    SaveCode = 1;
  end

  if ~isfield(File,'BestChains'),
    File.BestChains = zBestChains(File,1);
    SaveCode = 1;
  end

  if ~isfield(File,'Info') || ~isfield(File.Info,'ExpTechnique') || isempty(File.Info.ExpTechnique),
    File = zGetPDBInfo(File);          % get resolution and other info
    SaveCode = 1;
  end


  % obsolete code for reading .pdb1 files to get the biological unit coordinates; these can be obtained with symmetry operators in CIF files
  if 0 > 1,
    % ----------------- If it just read the coordinate file and no pairs, look at BUC

    if (ReadText == 1) && (isempty(strfind(PDBFilename,'.pdb1'))) && length(File.NT) > 0,
      bp = nnz(zSparseRange(abs(triu(File.Edge)),0,16));
                                                    % number of basepairs
      r  = bp / length(File.NT);                    % ratio of bp to nt
      if (r < 0.4) && isempty(strfind(lower(PDBFilename),'cif')),                                 % very few basepairs found
        if Verbose > 0,
          fprintf('Few basepairs found (%7.4f basepairs per nucleotide), reading the biological unit coordinates\n',r);
          zFlushOutput;
        end
        File1 = zGetNTData([PDBFilename '1'],4,1);   % read biological unit coords
        if length(File1.NT) > 0,
          bp = nnz(zSparseRange(abs(triu(File1.Edge)),0,16));
                                                    % number of basepairs
          r1 = bp / length(File1.NT);               % ratio of bp to nt
          if Verbose > 0,
            fprintf('Biological unit coordinates have %7.4f basepairs per nucleotide\n',r1);
            zFlushOutput;
          end
          if r1 > r,                              % pdb1 has higher basepair ratio
            File = File1;                         % use biological unit coords
          end
        end
      end
    end
  end

  File = orderfields(File);

  Saved = 0;

  if (ReadCode > 0) || (ClassifyCode > 0) || (SaveCode > 0),
%    zUpdatePDBInfo(File);
    zSaveNTData(File,Verbose);
    Saved = 1;
  end

  Files(f) = File;

end
