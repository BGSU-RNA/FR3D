% zSaveNTData saves File in File.Filename.analyzed
% Changed by Anton

function [] = zSaveNTData(Files,Verbose)

if nargin < 2,
  Verbose = 1;
end

for f=1:length(Files),

  File = Files(f);

  File.Modified = 0;                          % flag to know to save
  File.Distance = [];                         % clear; recompute on load

  if ~(exist('PrecomputedData','dir')),        % if directory doesn't yet exist
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

   [pathstr, name, ext] = fileparts(File.Filename); % Anton 5/30/2010
   destination = [pwd filesep 'PrecomputedData' filesep pathstr]; % Anton 5/30/2010
   if ~isempty(pathstr) && ~exist(destination,'dir') % Anton 5/30/2010
       mkdir(destination);
       path(path,destination);
   end   
   
   save([destination filesep name '.mat'],'File','-v6');  % Anton 5/30/2010, Anton 7/21/2011
%    save([pwd filesep 'PrecomputedData' filesep name '.mat'],'File'); %    Anton 5/30/2010
   if Verbose > 0,
     fprintf('Saved %s\n', [File.Filename '.mat']);
   end
  end

end
