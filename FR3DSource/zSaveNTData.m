% zSaveNTData(Files,Verbose) saves one or more files

function [] = zSaveNTData(Files,Verbose)

if nargin < 2,
  Verbose = 1;
end

for f=1:length(Files),
  File = Files(f);

  File.Modified = 0;                          % flag to know to save
  File.Distance = [];                         % clear; recompute on load

%  File.Filename = upper(File.Filename);      % probably dangerous to change case

  if File.NumNT >= 0,
    File.Pair = [];
    File.CI   = [];
    File.SizeCode = 1;

    for n=1:length(File.NT),
      File.NT(n).Loc = [];
    end

    [pathstr, name, ext] = fileparts(File.Filename); % Anton 5/30/2010
pathstr
name
ext

    if length(pathstr) == 0,                      % no path specified, so save in default location
      if ~(exist('PrecomputedData','dir')),        % if directory doesn't yet exist
        mkdir('PrecomputedData');
      end
      save([pwd filesep 'PrecomputedData' filesep upper(name) '.mat'],'File','-v6');  % Anton 5/30/2010, Anton 7/21/2011, Zirbel 2014/12/06 added upper()
      if Verbose > 0,
        fprintf('Saved %s\n', [File.Filename '.mat']);
      end
    else
      if ~exist(pathstr,'dir'),
        mkdir(pathstr);
        path(path,pathstr);
      end   
   
      save([pathstr filesep upper(name) '.mat'],'File','-v6');  % Anton 5/30/2010, Anton 7/21/2011, Zirbel 2014/12/06 added upper()
      if Verbose > 0,
        fprintf('Saved %s in %s\n', [name '.mat'], pathstr);
      end
    end
  end
end
