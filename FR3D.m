% FR3D conducts the search given in xSpecifyQuery in the PDB files listed
% in zFileNameList

% Change the list of PDB files to be searched by editing zFileNameList
% Change the query by editing xSpecifyQuery

if ~exist('GUIactive'),                      % if the GUI is not being used
  Query     = xSpecifyQuery;                 % get search parameters
  if isfield(Query,'SearchFiles'),           % if query specifies files
    Filenames = Query.SearchFiles;
  else
    Filenames = {'1s72'};                    % default
  end
end

% ----------------------------------------- Load PDB files if needed --------

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,2);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,2,File); % add PDB data if needed
end                       % SIndex tells which elements of File to search

% ------------------------------------------- Store actual filenames
%                                             rather than list name(s)

clear Filenames

for i = 1:length(SIndex),
  Filenames{i} = File(SIndex(i)).Filename;
end

% ------------------------------------------- Construct details of search ---

if isfield(Query,'Filename'),                % if query motif is from a file
  [File,QIndex] = zAddNTData(Query.Filename,0,File);  
                                             % load data for Query, if needed
  Query = xConstructQuery(Query,File(QIndex)); % preliminary calculations
else
  Query = xConstructQuery(Query);              % preliminary calculations
end

clear Search
Search.SaveName = [datestr(now,31) '-' Query.Name];  
                                  % use date and time to identify this search

%Query
%Search

if isfield(Query,'NumNT'),                    % if query is specified OK

% ------------------------------------------- Display query information------

fprintf('Query %s:', Query.Name);             % display query name

if isfield(Query,'Description'),
  fprintf(' %s\n', Query.Description);
else
  fprintf('\n');
end

% ------------------------------------------- Calc more distances if needed -

if Query.Geometric > 0,                       % if a geometric search
  tic                                         % keep track of time
  CalcFlag = 0;                               % if more distances were needed
  for f=1:length(SIndex),
    i = SIndex(f);
    if ceil(Query.DistCutoff) > ceil(max(max(File(i).Distance))),
      c = cat(1,File(i).NT(1:File(i).NumNT).Center);
      File(i).Distance = zMutualDistance(c,Query.DistCutoff); 
             % sparse matrix of center-center distances, up to Query.DistCutoff
      if length(File(i).NT) > 10,
%        zSaveNTData(File(i));                % don't bother, avoid mistakes
%        drawnow;
      end
      CalcFlag = 1;
    end
  end
  if CalcFlag > 0,
    fprintf('Calculated more distances in %5.3f seconds\n', toc);
  end
end

drawnow

% ------------------------------------------- Find candidates ---------------

starttime = cputime;

Candidates = xFindCandidates(File(SIndex),Query);  % screen for candidates

if ~isempty(Candidates),                         % some candidate(s) found
 if Query.Geometric > 0,
  [Discrepancy, Candidates] = xRankCandidates(File(SIndex),Query,Candidates);
  fprintf('Found %d candidates in the desired discrepancy range\n',length(Discrepancy));

   if (Query.ExcludeOverlap > 0) & (length(Discrepancy) > 0) ...
     & (Query.NumNT > 2),
     [Candidates, Discrepancy] = xReduceOverlap(Candidates,Discrepancy); 
                                                 % quick reduction in number
     [Candidates, Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,400); 
                                                % find top 400 distinct ones
     fprintf('Removed highly overlapping candidates, kept %d\n', length(Candidates(:,1)));
   end

 else
  A = [Candidates sum(Candidates')'];        % compute sum of indices
  N = Query.NumNT;                           % number of nucleotides
  [y,i] = sortrows(A,[N+1 N+2 1:N]);         % sort by this sum
  Candidates = Candidates(i,:);              % put all permutations together
  Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates
 end

% -------------------------------------------------- Save results of search

 Search.Query       = Query;
 Search.Filenames   = Filenames;
 Search.TotalTime   = cputime - starttime;
 Search.Date        = Search.SaveName(1:10);
 Search.Time        = Search.SaveName(12:18);
 Search.SaveName    = strrep(Search.SaveName,' ','_');
 Search.SaveName    = strrep(Search.SaveName,':','_');
 Search.Candidates  = Candidates;
 Search.Discrepancy = Discrepancy;

 Search = xAddFiletoSearch(File(SIndex),Search);

 if ~(exist('SearchSaveFiles') == 7),        % if directory doesn't yet exist
   mkdir('SearchSaveFiles');
 end

 save(['SearchSaveFiles' filesep Search.SaveName], 'Search');

% ------------------------------------------------ Display results

 fprintf('Entire search took %8.4f seconds, or %8.4f minutes\n', (cputime-starttime), (cputime-starttime)/60);

 if (~exist('GUIactive')) && (~isempty(Candidates)),
   xListCandidates(Search,Inf,1);
   Search = xDisplayCandidates(File(SIndex),Search);
   save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
 end

end

end