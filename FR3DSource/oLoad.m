% oLoad is a simple script which lists the saved searches in the SearchSaveFiles folder and prompts the user to load one.  This avoids the user having to copy and paste long filenames

more off                                   % turn off paging

d = dir([pwd filesep 'SearchSaveFiles']);

keep = ones(1,length(d));
for i = 1:length(d),
	if d(i).name(1) == '.',
		keep(i) = 0;
	end
	if ~isempty(strfind(d(i).name,'.txt')),
		keep(i) = 0;
	end
end
d = d(find(keep));

if length(d) > 0,

	if length(d) > 20,
		Filter = input('More than 20 files.  Enter filter text, or enter for none: ','s');
	else
		Filter = [];
	end

	keep = ones(1,length(d));
	for i = 1:length(d),
		if ~isempty(Filter),
			if isempty(strfind(lower(d(i).name),lower(Filter))),
				keep(i) = 0;
			end
		end
	end
	d = d(find(keep));

	if length(d) > 0,
		for i = 1:length(d),
			fprintf('%d   %s\n', i, d(i).name);
		end

		k = input('Please enter the number of the file you would like to load: ');

		load([pwd filesep 'SearchSaveFiles' filesep d(k).name]);
		Query = Search.Query;

		fprintf('This search found %d candidates\n',length(Search.Candidates(:,1)));

		oDisplay
	else
		fprintf('No files with filter text %s were found in SearchSaveFiles\n',Filter);
	end

else
	fprintf('No files found in SearchSaveFiles\n');
end
