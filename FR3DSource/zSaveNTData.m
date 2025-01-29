% zSaveNTData(Files,Verbose) saves one or more files

function [] = zSaveNTData(Files,Verbose)

if nargin < 2
  Verbose = 1;
end

for f=1:length(Files)
  File = Files(f);

  File.Modified = 0;                          % flag to know to save
  File.Distance = [];                         % clear; recompute on load

  %  File.Filename = upper(File.Filename);      % probably dangerous to change case

  if File.NumNT >= 0
    File.Pair = [];
    File.CI   = [];
    File.SizeCode = 1;

    for n=1:length(File.NT),
      File.NT(n).Loc = [];
    end

    [PDpath, name, ext] = fileparts(File.Filename);

    if length(PDpath) == 0,                      % no path specified, so save in default location
      PDpath = zPathToDirectory('PrecomputedData');
      if length(PDpath) == 0
        PDpath = zPathToDirectory('PDBFiles');
        if length(PFpath) == 0
          PDpath = pwd;
        end
      end
    end

    PDpath = [PDpath filesep 'PrecomputedData'];

    if ~exist(PDpath,'dir')
      mkdir(PDpath);
    end

    save([PDpath filesep upper(name) '.mat'],'File','-v6');  % Anton 5/30/2010, Anton 7/21/2011, Zirbel 2014/12/06 added upper()

    fprintf('Saved %s in %s\n', [File.Filename '.mat'], PDpath);

  end
end
