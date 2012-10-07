% xFR3DSearch conducts the search given in xSpecifyQuery in the PDB files 
% listed in zFileNameList

% Change the list of PDB files to be searched by editing zFileNameList
% Change the query by editing xSpecifyQuery

% The result of running xFR3DSearch is a variable Search with these fields:
%   
%    .Query                     % a description of the search parameters
%    .Filenames                 % names of files that were searched
%    .TotalTime                 % how much time the search took
%    .Date                      % date of the search
%    .Time                      % time the search was started
%    .SaveName                  % name for the saved search file
%    .Candidates                % L x N+1 matrix of indices of candidates
                                % L is the number of candidates found
                                % N is the number of search nucleotides
                                % the last column of Candidates is the index
                                % of the file in which the candidate was found
%    .Discrepancy               % geometric discrepancy of each candidate
                                % from the query, for geometric searches

if ~exist('Verbose'),
  Verbose = 1;                               % default is to print output
end

if ~exist('GUIactive') && ~exist('UsingLibrary'),  % FR3D should just search
  Query     = xSpecifyQuery;                 % get search parameters
  if isfield(Query,'SearchFiles'),           % if query specifies files
    Filenames = Query.SearchFiles;
  else
    Filenames = {'1s72'};                    % default
  end
end

% ----------------------------------------- Load PDB files if needed --------

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData(Filenames,0,[],Verbose);   % load PDB data
else
  [File,SIndex] = zAddNTData(Filenames,0,File,Verbose); %add PDB data if needed
end                       % SIndex tells which elements of File to search

% ------------------------------------------- Store actual filenames
%                                             rather than list name(s)

clear Filenames

for i = 1:length(SIndex),
  Filenames{i} = File(SIndex(i)).Filename;
end

% ------------------------------------------- Construct details of search ---

if ~exist('UsingLibrary'),
 if isfield(Query,'Filename'),                % if query motif is from a file
  [File,QIndex] = zAddNTData(Query.Filename,0,File);  
                                             % load data for Query, if needed
  Query = xConstructQuery(Query,File(QIndex)); % preliminary calculations
 else
  Query = xConstructQuery(Query);              % preliminary calculations
 end

end

clear Search
Search.SaveName = [datestr(now,31) '-' Query.Name];  
                                  % use date and time to identify this search

%Query
%Search

if isfield(Query,'NumNT'),                    % if query is specified OK

% ------------------------------------------- Display query information------

if Verbose > 0,
  fprintf('Query %s:', Query.Name);             % display query name

  if isfield(Query,'Description'),
    fprintf(' %s\n', Query.Description);
  else
    fprintf('\n');
  end
end

% ------------------------------------------- Calc more distances if needed -

for f=1:length(SIndex),
  i = SIndex(f);
  if isempty(File(i).Distance),
    dmin = 0;
  else
    dmin = ceil(max(max(File(i).Distance)));
  end

  if (ceil(Query.DistCutoff) > dmin) && (File(i).NumNT > 0),
    c = cat(1,File(i).NT(1:File(i).NumNT).Center);
    File(i).Distance = zMutualDistance(c,Query.DistCutoff); 
             % sparse matrix of center-center distances, up to Query.DistCutoff
  end
end

drawnow

% ------------------------------------------- Find candidates ---------------

starttime = cputime;

Candidates = xFindCandidates(File(SIndex),Query,Verbose);  % screen for candidates

if ~isempty(Candidates),                         % some candidate(s) found
 if Query.Geometric > 0,
  [Discrepancy, Candidates] = xRankCandidates(File(SIndex),Query,Candidates,Verbose);
  if Verbose > 0,
    fprintf('Found %d candidates in the desired discrepancy range\n',length(Discrepancy));
  end

   if (Query.ExcludeOverlap > 0) & (length(Discrepancy) > 0) ...
     & (Query.NumNT > 2),
     [Candidates, Discrepancy] = xReduceOverlap(Candidates,Discrepancy); 
                                                 % quick reduction in number
     [Candidates, Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,400); 
                                                % find top 400 distinct ones
     if Verbose > 0,
       fprintf('Removed highly overlapping candidates, kept %d\n', length(Candidates(:,1)));
     end
   end

 elseif Query.NumNT > 2,
  A = [Candidates sum(Candidates')'];        % compute sum of indices
  N = Query.NumNT;                           % number of nucleotides
  [y,i] = sortrows(A,[N+1 N+2 1:N]);         % sort by file, then this sum
  Candidates = Candidates(i,:);              % put all permutations together
  Discrepancy = (1:length(Candidates(:,1)))';% helps identify candidates
 else
  N = Query.NumNT;                           % number of nucleotides
  [y,i] = sortrows(Candidates,[N+1 1 2]);
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
  Search.SaveName    = strrep(Search.SaveName,'<','_');
  Search.SaveName    = strrep(Search.SaveName,'>','_');
  Search.SaveName    = strrep(Search.SaveName,'?','_');
  Search.SaveName    = strrep(Search.SaveName,'*','_');
  Search.SaveName    = strrep(Search.SaveName,'&','_');
  Search.Candidates  = Candidates;
  Search.Discrepancy = Discrepancy;

  Search = xAddFiletoSearch(File(SIndex),Search);

  if ~exist('UsingLibrary'),

    if ~(exist('SearchSaveFiles') == 7),     % if directory doesn't yet exist
      mkdir('SearchSaveFiles');
    end

    save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
 end
% ------------------------------------------------ Display results

 if Verbose > 0,
   fprintf('Entire search took %8.4f seconds, or %8.4f minutes\n', (cputime-starttime), (cputime-starttime)/60);
 end

 if (~exist('GUIactive')) && (~isempty(Candidates)) && ~exist('UsingLibrary'),
   xListCandidates(Search,Inf,1);
   Search = xDisplayCandidates(File(SIndex),Search);
   save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
 end

end

end