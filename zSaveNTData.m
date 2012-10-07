% zSaveNTData saves File in File.Filename.analyzed

function [void] = zSaveNTData(Files)

for f=1:length(Files),

  File = Files(f);

  File.Modified = 0;                          % flag to know to save
  File.Distance = [];                         % clear; recompute on load

  if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
    mkdir('PrecomputedData');
  end

  if File.NumNT >= 0,

   if File.SizeCode == 1,
     save([pwd filesep 'PrecomputedData' filesep File.Filename '.mat'],'File');
     fprintf('Saved %s\n', [File.Filename '.mat']);

     File = zSmallVersion(File);
   end

   save([pwd filesep 'PrecomputedData' filesep File.Filename '_small.mat'],'File');
   fprintf('Saved %s\n', [File.Filename '_small.mat']);
  end

end
