% oDisplay is a simple script that displays the current search.  It is useful in Octave, which does not have a GUI.

more off                                   % turn off paging

if exist('File') == 1 && exist('SIndex') == 1,
  Search = xDisplayCandidates(File(SIndex),Search);
else
	Search = xDisplayCandidates([],Search);
end
fprintf('Saving search\n')
zFlushOutput
save(['SearchSaveFiles' filesep Search.SaveName '.mat'], 'Search');
