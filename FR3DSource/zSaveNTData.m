% zSaveNTData saves File in File.Filename.analyzed

function [void] = zSaveNTData(Files,Verbose)

if nargin < 2,
  Verbose = 1;
end

for f=1:length(Files),

  File = Files(f);

  File.Modified = 0;                          % flag to know to save
  File.Distance = [];                         % clear; recompute on load

  if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
    mkdir('PrecomputedData');
  end

  File.Filename = upper(File.Filename);

  if File.NumNT >= 0,
   File.Pair = [];
   File.CI   = [];
   File.SizeCode = 1;

   for n=1:length(File.NT),
     File.NT(n).Loc = [];
   end

   save([pwd filesep 'PrecomputedData' filesep File.Filename '.mat'],'File');
   if Verbose > 0,
     fprintf('Saved %s\n', [File.Filename '.mat']);
   end
  end

end
