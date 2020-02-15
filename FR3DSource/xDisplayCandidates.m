% xDisplayCandidates(File,Search,Level) displays candidates
% graphically and allows various types of analysis

% It may be run directly from Matlab using the command:
%    [Search,File] = xDisplayCandidates([],Search);
% and after this,
%    [Search,File] = xDisplayCandidates(File,Search);

function [Search, File] = xDisplayCandidates(FullFile,Search,Level,UsingFull,Order,ShowNavWindow,Octave)

% if File is a text string (filename), load the file and display

if strcmp(class(FullFile),'char'),
  Filename = FullFile;
  FullFile = zAddNTData(Filename,0);
  File = FullFile;
end

if nargin < 7,
  Octave = 0;                                  % whether FR3D is being run through Octave
end

if exist('OCTAVE_VERSION') ~= 0,
  Octave = 1;
end

if ~isempty(strfind(path,'zirbel')),
%  Octave = 1;
end

if nargin < 6,
  ShowNavWindow = 0;                           % default is not to
end

if strcmp(class(Search),'char'),               % look at one or a few NTs
  S = zIndexLookup(File,Search);
  clear Search
  Search.Candidates = [S 1];
  Search.File = FullFile;
  UsingFull = 1;                             % use what was passed in
  File = FullFile;
  Search.Query.Geometric = 0;
  [L,N] = size(Search.Candidates);
  N = N - 1;                                   % number of nucleotides
  Search.Query.NumNT = N;
  Search.Discrepancy = 1:L;
  FIndex = 1:length(FullFile);
  for f = 1:length(File),
    Search.CandidateFilenames{f} = File(f).Filename;
  end
end

if strcmp(class(Search),'double'),
  S = Search;
  clear Search
  Search.Candidates = S;
  Search.File = FullFile;
  UsingFull = 1;                             % use what was passed in
  File = FullFile;
  Search.Query.Geometric = 0;
  [L,N] = size(Search.Candidates);
  N = N - 1;                                   % number of nucleotides
  Search.Query.NumNT = N;
  Search.Discrepancy = 1:L;
  FIndex = 1:length(FullFile);
  for f = 1:length(File),
    Search.CandidateFilenames{f} = File(f).Filename;
  end
end

if isempty(Search.Candidates)
  fprintf('There are no candidates to display\n');
  File = FullFile;
  return
end

[L,N] = size(Search.Candidates);
N = N - 1;                                   % number of nucleotides

if isfield(Search,'Limit'),
  Limit = Search.Limit;
else
  Limit = min(L,300);                          % for mutual discrep matrix
end

p = 1:L;                                     % default permutation for display
q(p) = 1:L;                                  % inverse permutation
                               % p : display position -> real candidate number
                               % q : real candidate number -> display position

if ~isfield(Search,'File'),
  UsingFull = 1;                             % use what was passed in
  File = FullFile;
  Search.Query.Geometric = 0;
  Search.Query.NumNT = N;
  Search.Discrepancy = 1:L;
  FIndex = 1:length(FullFile);
else
  UsingFull = 0;
  File = Search.File;                        % use what was saved w/ Search
  FIndex = 1:length(Search.File);
end

if ~isfield(Search,'SaveName'),
  Search.SaveName = 'these candidates';
end

if ~isfield(Search,'Query'),
  Search.Query.Geometric = 0;
  Search.Query.NumNT = length(Search.Candidates(1,:)) - 1;
end

Query = Search.Query;

if ~isfield(Search,'Discrepancy'),
  Search.Discrepancy = zeros(size(Search.Candidates(:,1)));
end

fontsize = 10;                               % for nucleotide numbers
NeighborhoodChanged = 0;

if nargin < 3,
	if Octave == 0,
	  MenuTitle  = 'Navigation options';
	else
		MenuTitle = ['Navigation options for ' Search.SaveName];
	end
  Level      = 0;
  QuitButton = 'Quit navigation';
else
  MenuTitle  = ['Subset depth ' num2str(Level)];
  QuitButton = 'Return to larger set';
end

if nargin < 5,
  if Query.Geometric == 0,
    Order = 2;
  else
    Order = 1;
  end
end

OrderText = {'by discrepancy from query', 'by file, then sum of nucleotide numbers', 'by similarity', 'by centrality', 'by pair criteria'};

if Octave > 0,
  DiaryFile = [Search.SaveName '.txt'];
else
  DiaryFile = 'diary.txt';
end

DiaryText = ['Append output to ' DiaryFile];
DiaryState = 0;

warning off

if ~isfield(Search,'Marked'),
  Search.Marked = zeros(1,L);         % allow marking certain candidates
end

if ~isfield(Search,'Disc'),
  Search.Disc         = sparse(zeros(L,L));
  Search.DiscComputed = sparse(zeros(1,L));
  if Query.Geometric > 0,
    Search.Disc(:,1) = Search.Discrepancy;
    Search.Disc(1,1) = 0;
    Search.DiscComputed(1,1) = 1;
  end
end

% ---------------- Make labels for candidates

for ii=1:L,
  f = Search.Candidates(ii,N+1);          % file number
  b = '';
  for j = 1:N,                            % list all bases
    b = [b File(f).NT(Search.Candidates(ii,j)).Base];
  end
  n = File(f).NT(Search.Candidates(ii,1)).Number;
  n = sprintf('%5s',n);
  if Search.Query.Geometric > 0,
      if isfield(Search,'AvgDisc'),
        d = sprintf('%6.4f',Search.AvgDisc(ii));
      else
        d = sprintf('%6.4f',Search.Discrepancy(ii));
      end
    else
      d = sprintf('%5d',Search.Discrepancy(ii)); % orig candidate number
    end
  Search.Lab{ii} = [File(f).Filename n ' ' b];

  if N == 2 && ~isempty(strfind(pwd,'zirbel')),
    ei = Search.Candidates(ii,1);
    ej = Search.Candidates(ii,2);
    et = strtrim(zEdgeText(File(f).Edge(ei,ej)));
    Search.Lab{ii} = [Search.Lab{ii} sprintf('%+5s',et)];
  end
end

% ----------------------------- find maximum gap between candidate nucleotides

[y,r] = sort(Search.Candidates(1,1:N)); % put nucleotides in increasing order

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(r,r),1);
else
  MaxDiff = Inf*ones(1,N-1);
end

maxinsert = zeros(1,N-1);
if N > 1,
  for c = 1:L,
    maxinsert = max(maxinsert,abs(diff(double(Search.Candidates(c,r))))-1);
  end
end

Display(1).p         = r;
Display(1).MaxDiff   = MaxDiff;
Display(1).MaxInsert = maxinsert;

if N == 2,
  Display(1).AtOrigin  = 1;           % put pairs in a standard orientation
else
  Display(1).AtOrigin = 0;
end

% ---------------------------- determine which conserved columns are on the same strand
% ---------------------------- number the columns by strand number

samestrand = eye(N,N);
for a = 1:N,
  for b = 1:N,
    if (isfield(Query,'Flank') && ~isempty(Query.Flank{a,b})) || (isfield(Query,'MaxDiffMat') && Query.MaxDiffMat(a,b) < Inf),
      samestrand(a,b) = 1;
      samestrand(b,a) = 1;
    end
    if max(abs(double(Search.Candidates(:,a))-double(Search.Candidates(:,b)))) <= 5,
      samestrand(a,b) = 1;
      samestrand(b,a) = 1;
    end
  end
end

samestrand = samestrand ^ N;

strandnumber = zeros(1,N);
for a = 1:N,
  if strandnumber(1,a) == 0,
    b = find(samestrand(a,:) > 0);
    strandnumber(1,b) = max(strandnumber) + 1;
  end
end

Display(1).strandnumber = strandnumber;

% --------- if there is no geometric model, align to the central candidate

if Query.Geometric == 0 || ~isfield(Query,'WeightedCenteredCenters'),
  [z,j]          = sort(sum(Search.Disc));           % sort by average discrepancy
  f              = Search.Candidates(j(1),N+1);
  Query.Indices  = double(Search.Candidates(j(1),1:N));
  Query.NT       = File(f).NT(Query.Indices);
  if ~isfield(Query,'LocWeight'),
    Query.LocWeight= ones(1,Query.NumNT);
  end
  Query          = xPrecomputeForDiscrepancy(Query);
  Query.Filename = 'Central candidate';
end

Display(1).rotation  = Query.NT(1).Rot;

% ------------------------------------------- Parameters to display candidates

Display = VisualizationOptions(Display,Octave,Search);    % set defaults

stop     = 0;                              % stop the menu?
i        = 1;                              % current window
nn       = 1;                              % current candidate

PlotMotif(File,Search,Query,Display,i);    % graph in display window i

if Octave == 0,
  rotate3d on
end
axis off
DisplayTable(File,Search,Query,Display,i)
drawnow

if N == 2,                                 % helpful for viewing basepairs
  view(2)
end

% ------------------------------- display menu -----------------------------

if ShowNavWindow > 0,
  if N == 2 && exist('xMutualIDI') ==2,              % 2-NT candidates
    figure(98)
    axis([1 Limit+1 1 Limit+1]);
  end

  figure(99)
  axis([1 Limit+1 1 Limit+1]);
end

while stop == 0,

 if ShowNavWindow > 0,
  % ---------------------------------------- Display table of discrepancies
  figure(99)
  ax = axis;
  clf
  pp = p(1:Limit);

  if Octave == 1,
	  % --------------------------- Mark displayed candidates with blue on diagonal
  	for ww = 1:length(pp),
	  	Search.Disc(pp(ww),pp(ww)) = 0;
	  end
		ww = q(find(Search.Marked));
	  for j = 1:length(ww),
	  	Search.Disc(pp(ww(j)),pp(ww(j))) = 0.4;
	 	end
	  for j = 1:length(Display),
	  	Search.Disc(pp(q(Display(j).n)),pp(q(Display(j).n))) = 1;
	  end
	  zGraphDistanceMatrix(Search.Disc(pp,pp),Search.Lab(pp));
	else
	  % --------------------------- Mark displayed candidates with white symbol
	  zGraphDistanceMatrix(Search.Disc(pp,pp),Search.Lab(pp));
	  hold on
	  co = {'w*','wo','wd','ws','wv','w<','w>','w^','w+','wx'};
	  co = [co co co co co co co co];
	  if ~exist('nowhitemark.txt'),
	    for j = 1:length(Display),
	      plot(q(Display(j).n)+0.5,q(Display(j).n)+0.5,co{j});
	    end
	    m = q(find(Search.Marked));
	    plot(m+0.5,m+0.5,'w.');
	  end
	end

  if Display(1).SimilarityToggle == 1,
%    SU = Display(1).SimilarityUnique + 1;
%    hold on
%    plot([1 Limit],[SU SU],'w');
%    plot([SU SU],[1 Limit],'w');
  end


  if Limit < L,
    title(['Discrepancies between first ' num2str(Limit) ' candidates, ordered by ' OrderText{Order}]);
  else
    title(['Discrepancies between all candidates, ordered by ' OrderText{Order}]);
  end

  colormap('default');
  map = colormap;
  map = map((end-8):-1:8,:);
  colormap(map);
  caxis([0 1.4]);

  if Octave == 0,
	  colorbar('location','eastoutside');
    set(gcf,'Name','Navigation window; click here, then click the "Navigate" button');
  else
    set(gcf,'Name','Navigation window; select Navigate with Figure 99 then click here');
	end

  if N == 2 && exist('xMutualIDI') == 2,              % 2-NT candidates

	  figure(98)
	  ax = axis;
	  clf
	  pp = p(1:Limit);
	  zGraphDistanceMatrix(Search.IDI(pp,pp));
	%  zGraphDistanceMatrix(Search.IDI(pp,pp),Search.Lab(pp));
	  hold on
	  co = {'w*','wo','wd','ws','wv','w<','w>','w^','w+','wx'};
	  co = [co co co co co co co co];
	  for j = 1:length(Display),
	    plot(q(Display(j).n)+0.5,q(Display(j).n)+0.5,co{j});
	  end
	  m = q(find(Search.Marked));
	  plot(m+0.5,m+0.5,'w.');
	%  axis(ax);
	  if Limit < L,
	    title(['IsoDiscrepancies between first ' num2str(Limit) ' candidates, ordered by ' OrderText{Order}]);
	  else
	    title(['IsoDiscrepancies between all candidates, ordered by ' OrderText{Order}]);
	  end
	  colormap('default');
	  map = colormap;
	  map = map((end-8):-1:8,:);
	  colormap(map);
	  caxis([0 5]);
	  if Octave == 0,
		  colorbar('location','eastoutside');
		end
	  set(gcf,'Name','Navigation window; click here, then click the "Navigate" button');

	  fprintf('Counts of base combinations found in this set.\n');

	  counts = zeros(4,4);
	  for i = 1:L,
	    f = Search.Candidates(i,3);
	    a = Search.Candidates(i,1);
	    b = Search.Candidates(i,2);
	    c1 = File(f).NT(a).Code;
	    c2 = File(f).NT(b).Code;
	    counts(c1,c2) = counts(c1,c2) + 1;
	  end

	  Letters = 'ACGU';

	  fprintf('        A      C      G      U\n');
	  for i = 1:4,
	    fprintf('%c   %5d  %5d  %5d  %5d\n', Letters(i), counts(i,1), counts(i,2), counts(i,3), counts(i,4));
	  end
	  fprintf('\n');

  end
 end

  Buttons = {'Next candidate','Previous candidate', ...       % 1,2
         'Add plot', DiaryText, ...                           % 3,4
         'Visualization options','Jump to candidate', ...     % 5,6
         'Mark/Unmark current','Reverse all marks', ...       % 7,8
         'Display marked only', ...                           % 9
         'List to screen','Write to PDB', ...                 % 10,11
         'Sort by centrality', 'Order by Similarity', ...     % 12,13
         'Show Alignment', ...                                % 14
         'Show Scatterplot', 'Navigate with Figure 99', ...   % 15, 16
         'Rotate 20 degrees',QuitButton};                     % 17,18

  if Octave > 0,
    k = oMenu(MenuTitle,Buttons);
  else
    k=menu(MenuTitle,Buttons);
  end

  currfig = gcf;
  try
    ii = currfig.Number;
  catch
      ii=gcf;                                 % get current active figure
  end
  if (abs(ii) > length(Display)) || (ii == 0), % other window active?
    ii = i;
  end
  i = ii;                                 % record and save active figure
  i = min(i,length(Display));

  figure(i)
  [az,el]=view;                          % get current view (orientation)
  Display(1).az = az;
  Display(1).el = el;
  Display(1).x=xlim;                     % current x, y, z limits
  Display(1).y=ylim;
  Display(1).z=zlim;

  % ------------------------------------------- Want the navigation window?

  if any(k == [12 13 16]),
    ShowNavWindow = min(2,1+ShowNavWindow);
  end

  if k ~= 12,
    Display(1).Centrality = 0;              % user has not just selected this
  end

 % ------------------------------------------- Calculate distance matrix

 if ShowNavWindow == 1,                     % just indicated to show this
   if Limit < L,
     fprintf('Calculating discrepancies between first %d candidates\n',Limit);
   else
     fprintf('Calculating discrepancies between all %d candidates\n',L);
   end
   drawnow
   Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies

   % ------------------------------------------- Calculate IDI matrix

   if N == 2 && exist('xMutualIDI') == 2,              % 2-NT candidates

     Search = xMutualIDI(File,Search,Limit); % calculate some discrepancies

     for ii=1:L,
       f = Search.Candidates(ii,N+1);          % file number
       b = '';
       for j = 1:min(4,N),
         b = [b File(f).NT(Search.Candidates(ii,j)).Base];
       end
       n = File(f).NT(Search.Candidates(ii,1)).Number;
       n = sprintf('%4s',n);
       if Search.Query.Geometric > 0,
         if isfield(Search,'AvgDisc'),
           d = sprintf('%6.4f',Search.AvgDisc(ii));
         else
           d = sprintf('%6.4f',Search.Discrepancy(ii));
         end
       else
         d = sprintf('%5d',Search.Discrepancy(ii)); % orig candidate number
       end
       Search.Lab{ii} = [b n ' ' File(f).Filename];
     end
   end
 end

  % ------------------------------------------- respond to menu choice

  switch k                                      % k is the menu choice
    case 1                                      % next plot
      n = Display(i).n;                         % actual candidate displayed

      if (ShowNavWindow > 0) && (q(n) + 1 > Limit) && (Limit < L),
        Limit = min(Limit*2,L);                 % increase limit
        Search.Limit = Limit;
        fprintf('Increased display limit to %d; calculating more discrepancies\n',Limit);
        drawnow
        Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies

        if ~isempty(strfind(pwd,'zirbel')),
          fprintf('Ordering by similarity\n');
          tic
%          p(1:Limit) = TSPGreedyInsertion(Search.Disc(1:Limit,1:Limit),[],100);
          p(1:Limit) = OLO(Search.Disc(1:Limit,1:Limit),'average');
          toc
        else
          p(1:Limit) = zOrderbySimilarity(Search.Disc(1:Limit,1:Limit));
        end

        p((Limit+1):L) = (Limit+1):L;
        q(p) = 1:L;                             % inverse permutation
      elseif q(n) + 1 > L,                      % q(n) is display order
        Display(i).n = p(1);
      else
        Display(i).n = p(q(n) + 1);
      end

    case 2                                      % Previous Plot
      n = Display(i).n;                         % actual candidate displayed
      if q(n) - 1 < 1,                          % q(n) is display order
        Display(i).n = p(L);
      else
        Display(i).n = p(q(n) - 1);
      end

    case 3                                      % Add plot
      Display(end+1) = Display(i);              % use current settings
      i = length(Display);                      % current figure number
      figure(i);

    case 4
      if DiaryState == 0,
        diary(['SearchSaveFiles' filesep DiaryFile]);
        DiaryState = 1;
        DiaryText = 'Stop appending output';
      else
        diary off
        DiaryState = 0;
        DiaryText = ['Append output to ' DiaryFile];
      end

    case 5                                      % toggle sugar & superimpose
      [Display, NeighborhoodChanged] = VisualizationOptions(Display,Octave);

    case 6
      jtc = input('Enter the PDB ID or nucleotide number or base string: ','s');
      jtc = strrep(jtc,'''','');
      keep = zeros(L,1);
      for jt = 1:L,
        if ~isempty(strfind(lower(Search.Lab{jt}),lower(jtc))),
          keep(jt,1) = 1;
        end
      end
      if sum(keep) == 1,
        Display(i).n = find(keep);
      elseif sum(keep) == 0,
        fprintf('No matching candidate found\n');
      else
        fprintf('Multiple matches:\n');
        fk = find(keep);
        for jt = 1:length(fk),
          fprintf('%d) %s\n',jt,Search.Lab{fk(jt)});
        end
        try
          jt = input('Enter number to jump to it: ','s');
          jt = strrep(jt,'''','');
          jt = str2num(jt);
          if ~isempty(jt) && jt > 0 && jt <= length(fk),
            Display(i).n = fk(jt);
          end
        end
      end

    case 7                                      % mark/unmark current cand
      Search.Marked(Display(i).n) = 1-Search.Marked(Display(i).n); %toggle
      n = Display(i).n;                         % actual candidate displayed
      if q(n) + 1 > L,                          % q(n) is display order
        Display(i).n = p(1);
      else
        Display(i).n = p(q(n) + 1);
      end

    case 8                                      % reverse all marks
      Search.Marked = 1-Search.Marked;

    case 9                                      % display marked only
      j = find(Search.Marked);
      if length(j) > 0,
        [y,m] = sort(q(j));
        j = j(m);                               % put j in display order
        Search2 = SearchSubset(Search,j);
        xDisplayCandidates(File(FIndex),Search2,Level+1,UsingFull,Order,ShowNavWindow);
      end

    case 10                                      % list on screen
      j  = find(Search.Marked);
      jj = find(Search.Marked == 0);
      if (length(j) > 0) && (length(jj) > 0),
        [y,m] = sort(q(j));
        j = j(m);                                % put j in display order
        Search2 = SearchSubset(Search,j);
        fprintf('Marked candidates listed first\n');
        drawnow
        xListCandidates(Search2,Inf);

        [y,m] = sort(q(jj));
        jj = jj(m);                             % put jj in display order
        Search2 = SearchSubset(Search,jj);
        fprintf('Unmarked candidates listed second\n');
        drawnow
        xListCandidates(Search2,Inf);
      else
        Search2 = SearchSubset(Search,p);
        xListCandidates(Search2,Inf);
      end

    case 11                                     % write PDB of all
      SearchT = Search;
      Search = SearchSubset(Search,p);
      xWriteCandidatePDB(Search);
      if Level > 0,
        SN = [Search.SaveName '_Subset_' datestr(now,31)];
        SN    = strrep(SN,' ','_');
        SN    = strrep(SN,':','_');
        Search.SaveName = SN;
        save(['SearchSaveFiles' filesep SN], 'Search');
      end
      Search = SearchT;

      if exist('xListCandidatesByID'),
        xListCandidatesByID(Search);
      end

    case 12                                     % sort by centrality

      if (Display(1).Centrality == 1) && (min(Limit*2,L) > Limit),
        Limit = min(Limit*2,L);                 % increase limit
        Search.Limit = Limit;
        fprintf('Increased display limit to %d; calculating more discrepancies\n',Limit);
        drawnow
        Search = xMutualDiscrepancy(File,Search,Limit); % calculate some discrepancies
      end

      [z,j] = sort(max(Search.Disc(1:Limit,1:Limit)));% sort by max discrepancy
%      [z,j] = sort(sum(Search.Disc));           % sort by average discrepancy
      S.AvgDisc  = z / (Limit - 1);             % average discrep among these
      p(1:Limit) = j;
      p((Limit+1):L) = (Limit+1):L;
      q(p) = 1:L;
%      Search = xSortByCentrality(File(FIndex),Search,Level,UsingFull);
      Order = 4;
      Display(1).Centrality = 1;                % user has just clicked this

    case 13                                     % order by similarity
      Display(1).SimilarityToggle = 1 - Display(1).SimilarityToggle;

      if Limit > 1,
      	if Octave == 1,
      		fprintf('Ordering by similarity\n');
          zFlushOutput;
      	end

        if ~isempty(strfind(pwd,'zirbel')),
%          p(1:Limit) = TSPGreedyInsertion(Search.Disc(1:Limit,1:Limit),[],100);
          p(1:Limit) = OLO(Search.Disc(1:Limit,1:Limit),'average');
        else
          p(1:Limit) = zOrderbySimilarity(Search.Disc(1:Limit,1:Limit));
        end

        p((Limit+1):L) = (Limit+1):L;
      else
        p = 1;
      end

      q(p) = 1:L;
      Order = 3;

      if N == 2 && exist('xMutualIDI') == 2,              % 2-NT candidates
        p(1:Limit) = zOrderbySimilarity(Search.IDI(1:Limit,1:Limit));
        p((Limit+1):L) = (Limit+1):L;
        q(p) = 1:L;
        Order = 3;
      end

      if Display(1).SimilarityToggle > 0,

        % go through all candidates, see if there are two of each, or three of each, if so, double or triple the limit, then order by similarity, then go through and keep just the first instance of each distinct candidate.

        SeenOnce = zeros(1,length(p));          %
        SortedCand = Search.Candidates(p,:);
        SortedCand(:,1:N) = sort(SortedCand(:,1:N),2);  % sort indices
        CandDist = zDistance(double(SortedCand));

%        figure(32)
%        clf
%        spy(CandDist == 0);

        fprintf('Displaying the first instance of each candidate\n');
        zFlushOutput;

        clear order

        for m = 1:Limit,
          if any(CandDist(m,1:(m-1)) == 0),      % not the first one
            order(m) = m + 100000;
          else
            order(m) = m;
          end
        end
        for m = (Limit+1):length(p),
          order(m) = m + 200000;
        end

        [o,m] = sort(order);

        q = p(m);                               % re-order

        p = q;                                  % new ordering

        m = length(find(order < 100000));       % number of distinct instances

        fprintf('Found %d distinct instances\n',m);
				zFlushOutput;

        q = zOrderbySimilarity(Search.Disc(p(1:m),p(1:m)));

        r = p(1:m);
        r = r(q);

        p(1:m) = r;                             % re-order distinct instances

        q(p) = 1:L;

        Display(1).SimilarityUnique = m;

        Search.Marked(p(1:m)) = 1;

        if length(j) > 0,
          j = p(1:m);
          Search2 = SearchSubset(Search,j);
          xDisplayCandidates(File(FIndex),Search2,Level+1,UsingFull,Order,ShowNavWindow);
        end
      end

    case 14                                     % align
      Text = xAlignCandidates(File(FIndex),Search,1,p);
%      Text = xFASTACandidates(File(FIndex),Search,1);
      for t = 1:length(Text),
        fprintf('%s\n', Text{t});
      end

    case 15                                     % scatterplot
      ViewParam.Color  = 6;
      ViewParam.FigNum = length(Display)+1;
      ViewParam.Normal = 0;
      ViewParam.ClassLimits = 1;
      ViewParam.Order = p;
      ViewParam.Color = 2;                      % color by interaction; formerly 1, color by position in list
      ViewParam.Octave = Octave;
      p = xScatterPairs(Search,1,2,ViewParam);
      q(p) = 1:L;
      Order = 5;

    case 16
      if ShowNavWindow == 2,                    % already displayed
      	if Octave == 1,
          fprintf('Click near the diagonal to select a candidate, away from the diagonal to mark a range of candidates\n');
          figure(99)
      		[x,y,b] = ginput(1);
      		pt = [x y];
      	else
          figure(99)
	        pt = get(gca,'CurrentPoint');
	      end
      else
        pt(1,1) = Display(i).n;                 % current candidate
        pt(1,2) = Display(i).n;
      end

      if abs(pt(1,1)-pt(1,2)) > Limit/20,           % clicked off the diagonal
        Search.Marked = 0 * Search.Marked;      % unmark all candidates
        a = sort(pt(1,[1 2]));
        j = p(max(1,floor(a(1))):min(L,floor(a(2))));      %
        Search.Marked(j) = ones(1,length(j));   % select these candidates
      else                                      % clicked near the diagonal
        newn = max(min(floor(pt(1,1)),L),1);
        Display(i).n = p(newn);
      end

    case 17
      deg = -20*pi/180;
      Display(1).rotation = Display(1).rotation * [cos(deg) 0 sin(deg); 0 1 0; (-sin(deg)) 0 cos(deg)];

    case 18                                     % quit Display
      if exist('fidOUT','var')
        fclose(fidOUT);
      end
      diary off
      stop = 1;

  end  % switch statement for menu

  if k ~= 13,
    Display(1).SimilarityToggle = 1;
  end

  if any([1 2 3 7 16 17] == k),
      PlotMotif(File(FIndex),Search,Query,Display,i);
     	axis off
  end

  if (Display(i).n ~= nn) || (k == 4) || (k == 5 && NeighborhoodChanged) || (L == 1 && k <= 2),
    DisplayTable(File(FIndex),Search,Query,Display,i)
    nn = Display(i).n;
    NeighborhoodChanged = 0;
  end

  if any([4 5 6 8] == k),
    for j=1:length(Display)
      PlotMotif(File(FIndex),Search,Query,Display,j);
      axis off
    end
  end

  if length(Display) > 1,
    zLinkFigures(1:length(Display));
  end

  if ShowNavWindow == 1,
    ShowNavWindow = 2;
  end

  figure(i)
  if Octave == 0,
    rotate3d on
  end
  drawnow

end  % end while

if UsingFull == 0,
  File = FullFile;
  FullFile = [];
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function PlotMotif(File,Search,Query,Display,i)

  N = Query.NumNT;

  figure(i)
  clf

  if (Display(1).superimpose == 1),
    Indices = Query.Indices;
    if (Query.NumNT > 2),
      R = eye(3);
      S = mean(cat(1,Query.NT.Center));
    else
      R = Query.NT(1).Rot;
      S = mean(cat(1,Query.NT.Center));
    end
    MVP.Rotation      = R * Display(1).rotation;
    MVP.Shift         = S;
    MVP.LineStyle     = '-.';
    MVP.LineThickness = 1;
    MVP.Sugar         = Display(1).superbackbone;
    MVP.ConnectSugar  = 0;
    MVP.Grid          = 0;
    MVP.LabelBases    = Display(1).labelbases;
    zDisplayNT(Query,1:N,MVP);
  end

  [s,t] = size(Search.Candidates);
  n       = Display(i).n;
  f       = Search.Candidates(n,N+1);
  Indices = double(Search.Candidates(n,1:N));

  if isfield(File(f),'Filename'),
    FN = File(f).Filename;
  else
    FN = '';
  end

  nt = File(f).NT(Indices(1));
  Title = [strrep(FN,'_','\_') ' ' nt.Base nt.Number];
  for j=2:min(10,length(Indices)),
    nt = File(f).NT(Indices(j));
    Title = [Title '-' nt.Base nt.Number];
  end;

  VP.Sugar    = Display(1).displaybackbone;
  VP.LabelBases = Display(1).labelbases;
  VP.ShowBeta   = Display(1).showbeta;

  if Query.NumNT > 2,                                  % align to the model
    MC = Query.CenteredCenters;                        % nucleotide centers minus the overall center
    CandiCenters = cat(1,File(f).NT(Indices).Center);  % centers for this candidate
    CC = CandiCenters - ones(N,1)*Query.LocWeight * CandiCenters / sum(Query.LocWeight);  % centered centers for this candidate

    R = zBestRotation(MC, CC, Query.LocWeight);
    S = mean(CandiCenters);
  elseif length(Indices) > 1,
    R = File(f).NT(Indices(1)).Rot;
    S = mean(cat(1,File(f).NT(Indices).Center));
  else
    R = File(f).NT(Indices(1)).Rot;
    S = File(f).NT(Indices(1)).Center;
  end

  VP.Rotation = R * Display(1).rotation;
  VP.Shift    = S;
  VP.Grid     = 0;

  if exist('amal.txt','file') > 0 && N == 3,
     VP.AtOrigin = 2;
  end

  if ~isempty(strfind(pwd,'zirbel')) && N == 2 && Display(1).AtOrigin == 1,
     VP.AtOrigin = 1;
  end

  OrigIndices = Indices;

  if max(Display(1).neighborhood) > 0,
    NeighIndices = xNeighborhood(File(f),Indices,Display(1).neighborhood,Display(1).strandnumber);
    Indices = [Indices NeighIndices];
  else
    NeighIndices = [];
  end

%  if Display(1).colorstate == 0 && max(Display(1).neighborhood) > 0,
%    VP.NucleotideBrightness(1:length(Indices)) = 1;
%    VP.NucleotideBrightness((length(Indices)+1):(length(Indices)+length(NeighIndices))) = 0.5;
%  end

  strandnumber = [Display(1).strandnumber zeros(1,length(NeighIndices))];
  NN = length(strandnumber);

  if Display(1).colorstate == 1,

    samestrand = eye(NN,NN);
    for a = 1:NN,
      for b = (a+1):NN,
        if abs(double(Indices(a)) - double(Indices(b))) <= 3,
          samestrand(a,b) = 1;
          samestrand(b,a) = 1;
        end
      end
    end

    samestrand = samestrand ^ NN;

    % sort indices by strand number and then by distance from 1st nucleotide?  that might color them more consistently

    for a = 1:NN,
      b = find(samestrand(a,:) > 0);     % find others in the same strand
      if strandnumber(1,a) > 0,            % already assigned to a strand
        strandnumber(1,b) = strandnumber(1,a); % assign to the same strand as a
      else
        strandnumber(1,b) = max(strandnumber) + 1;
      end
    end

%    [y,i] = sort(Indices);
%    Indices = Indices(i);
%    strandnumber = strandnumber(i);

    Palette = [255 0 51; ...     % red
               205 173 0; ...    % gold
               51 153 255; ...   % lt blue
               153 30 255; ...   % purple
               139 69 19; ...    % saddlebrown
               51 0 255; ...     % bk blue
               51 204 51; ...    % green
               128 128 128; ...  % gray
               255 102 0; ...    % orange
               ];

    if max(strandnumber) > length(Palette(:,1)),
      fprintf('Too many strands for the colors that have been defined, re-using some.\n');
    end

    % color the strands from light to dark

    for j = unique(strandnumber),
      a = find(strandnumber == j);

      [y,i] = sort(Indices(a));
      a = a(i);

      jj = 1 + mod(j-1,length(Palette(:,1)));       % re-use palette colors if necessary

      p = Palette(jj,:)/255;
      mi = 0.6;                                     % minimum multiplier
      ma = 1/max(p);                              % maximum multiplier

      if length(a) > 1,
        step = (mi-ma)/(length(a)-1);
        step = max(-0.2,step);                      % no need for drastic brightness changes
        for b = 1:length(a),
          VP.Colors(a(b),:) = p * (ma + step*(b-1));
          VP.NumberColors(a(b),:) = [0 0 0];
        end
      else
        VP.Colors(a,:) = p;
        VP.NumberColors(a,:) = [0 0 0];
      end
    end

    VP.Colors = min(1,VP.Colors);
    VP.Colors = max(0,VP.Colors);
  end

  VP.LineThicknesses = 2*ones(1,NN);
  if max(Display(1).neighborhood) > 0,
    VP.LineThicknesses(1:N) = 5;
  end

  VP.BackboneTrace = Display(1).backbonetrace;
  VP.AABackboneTrace = Display(1).aabackbonetrace;

  zDisplayNT(File(f),Indices,VP);

%  xlabel([File(f).Filename ' ' File(f).Info.Descriptor]);
  fprintf([File(f).Filename ' ' File(f).Info.Descriptor '\n']);

  if Display(1).nearbyatoms > 0,
    if isfield(File(f),'AA') && length(File(f).AA) > 0,
      c = cat(1,File(f).NT(OrigIndices).Center);      % nucleotide centers
      clear a
      for b = 1:length(File(f).AA),
        if ~isempty(File(f).AA(b).Loc),
          a(b,:) = File(f).AA(b).Center;               % amino acid centers
        else
          a(b,:) = [Inf Inf Inf];
        end
      end

      D = zDistance(c,a);
      [i,j,k] = find(D);                          %

      aadist = max(Display(1).neighborhood(8),8);      % within 8 Angstroms
      w = find(k <= aadist);

  %  	try
      if ~isempty(w),
        zDisplayAA(File(f),unique(j(w)),VP);
      end
  %   end
    else
      fprintf('No amino acids in structure %s\n',File(f).Filename);
    end
  end




  set(gcf,'Name',strrep(Title,'\_','_'));

  if isfield(Search,'AvgDisc'),
    xlabel(['Candidate ',int2str(n),' of ',int2str(s),'   Average discrepancy from others ', num2str(Search.AvgDisc(n))]);
  elseif Query.Geometric > 0,
    xlabel(['Candidate ',int2str(n),' of ',int2str(s),'   Discrepancy ',...
          num2str(Search.Discrepancy(n))]);
  else
    xlabel(['Candidate ',int2str(n),' of ',int2str(s)]);
  end

  if Search.Marked(n) == 1;
    yl = 'Marked';
  else
    yl = '';
  end
  ylabel(yl);

  axis equal
  axis vis3d
%  view(2)
  view([Display(1).az Display(1).el]);
  drawnow
%  axis off

  if exist('amal.txt','file') > 0 && N == 3 && Display(1).neighborhood == 0,
    zAmalTripleStudy    % Commands for Amal's study of base triples:
  end

%    set(gcf,'Renderer','OpenGL')
%    set(gcf,'Renderer','zbuffer')

% ------------------------------------------------- Display table

function  DisplayTable(File,Search,Query,Display,i)

    N       = Query.NumNT;
    n       = Display(i).n;
    f       = Search.Candidates(n,N+1);
    Indices = double(Search.Candidates(n,1:N));

    if isfield(Search,'AvgDisc'),
      fprintf('Average discrepancy from others %6.4f', Search.AvgDisc(n));
    elseif Query.Geometric > 0,
      fprintf('Discrepancy %6.4f', Search.Discrepancy(n));
    else
      fprintf('Candidate #%d of %d', Search.Discrepancy(n), length(Search.Candidates(:,1)));  % integer is cand num
    end

    if Search.Marked(n) == 1;
      fprintf('Marked\n');
    end

    if max(Display(1).neighborhood) > 0,
      Indices = sort([Indices xNeighborhood(File(f),Indices,Display(1).neighborhood,Display(1).strandnumber)]);
    end

    zShowInteractionTable(File(f),double(Indices));

    if isfield(File(f),'BasePhosphate'),
      zBasePhosphateTable(File(f),double(Indices));
    end
    if isfield(File(f),'BaseRibose'),
      zBaseRiboseTable(File(f),double(Indices));
    end

    drawnow


% -------------------------------------------------- Select subset of search

function [Search2] = SearchSubset(Search,j)

  Search2             = Search;
  Search2.Candidates  = Search.Candidates(j,:);
  Search2.Discrepancy = Search.Discrepancy(j);
  Search2.Marked      = Search.Marked(j);
  Search2.Disc        = Search.Disc(j,j);
  Search2.DiscComputed= Search.DiscComputed(1,j);
  Search2.Limit       = length(j);
  if isfield(Search,'Lab'),
    Search2.Lab         = Search.Lab(j);
  end

% -------------------------------------------------- Set visualization options
function [Display,NeighborhoodChanged] = VisualizationOptions(Display,Octave,Search)

NeighborhoodChanged = 0;

fontsize = 10;

if ~isfield(Display,'colorstate'), % set default values
  Display(1).n                = 1;     % which candidate is in display window 1
  Display(1).colorstate       = 0;     % 0 will color A red, etc., 1 color by strand
  Display(1).displaybackbone  = 1;     % display sugars or not
  Display(1).neighborhood     = zeros(1,8);     % what neighborhood to show
  Display(1).superimposestate = 0;     % superimpose the first candidate? backbone too?
  Display(1).superimpose      = 0;     % superimpose the first candidate?
  Display(1).superbackbone    = 0;     % show sugars of first when superimposing?
  Display(1).labelbases       = 10;    % show nucleotide numbers
  Display(1).az               = -37.5; % standard view
  Display(1).el               = 30;
  Display(1).az               = 0; % standard view
  Display(1).el               = 0;
  Display(1).toggle           = 1;     % default view, see below
  Display(1).nearbyatoms      = 0;     % default is to not show waters, amino acids
  Display(1).showbeta         = 0;     % default is to not show beta factors
  Display(1).backbonesuper    = 1;     % superimposing and sugars/backbones
  Display(1).SimilarityToggle = 1;     % 1 means eliminate duplicates; 0 not
  Display(1).SimilarityUnique = 1;
  Display(1).Centrality       = 0;     % user has not just clicked sort by centrality
  Display(1).backbonetrace    = 0;     % don't draw a separate backbone trace
  Display(1).aabackbonetrace  = 0;     % don't draw a separate backbone trace
  if isfield(Search,'oExploreNT'),
    Display(1).colorstate = 1;
%    Display(1).neighborhood = [1 0 0 1 1 0 0 12];
    Display(1).backbonetrace = 1;        % don't draw a separate backbone trace
  end
else

  k = 0;

% while k ~= -1,

    switch Display(1).colorstate,
    case 0,
      NextColorState = 'Color by strand';
    case 1,
      NextColorState = 'Color A red, C yellow, G green, U blue';
    end

    switch Display(1).displaybackbone,
    case 0,
      NextBackboneState = 'Display backbone';
    case 1,
      NextBackboneState = 'Do not display backbone';
    end

    switch Display(1).labelbases,
    case 0,
      NextLabelBasesState = 'Label bases';
    otherwise,
      NextLabelBasesState = 'Do not label bases';
    end

    switch Display(1).superimposestate,
    case 0,
      NextSuperimposeState = 'Superimpose query and backbone';   % 1
    case 1,
      NextSuperimposeState = 'Superimpose query, no backbone';   % 2
    case 2,
      NextSuperimposeState = 'Do not superimpose query';         % 0
    end

    switch Display(1).nearbyatoms,
    case 0,
      NextNearbyAtomState = 'Show amino acids';
    case 1,
      NextNearbyAtomState = 'Do not show amino acids';
    end

    Buttons = {...
              NextColorState, ...                                % 1
              NextBackboneState, ...                             % 2
              NextLabelBasesState, ...                           % 3
              NextSuperimposeState, ...                          % 4
              NextNearbyAtomState, ...                           % 5
              'Neighborhood: none', ...                          %  6, all zeros
              'Neighborhood: fill in strands', ...               %  7, 1 in position 1
              'Neighborhood: extend strands 1 nucleotide', ...   %  8, 1 in position 2
              'Neighborhood: extend strands 2 nucleotides', ...  %  9, 1 in position 3
              'Neighborhood: extend strands 3 nucleotides', ...  % 10, 1 in position 4
              'Neighborhood: basepairs made', ...                % 11, 1 in position 5
              'Neighborhood: stacks made', ...                   % 12, 1 in position 6
              'Neighborhood: other interactions made', ...       % 13, 1 in position 7
              'Neighborhood: 8A', ...                            % 14, 8 in position 8
              'Neighborhood: 12A', ...                           % 15, 12 in position 8
              'Neighborhood: 16A', ...                           % 16, 16 in position 8
              'Backbone trace', ...                              % 17 trace the backbone
              };

    MenuTitle = 'Display options';

    if Octave > 0,
      k = oMenu(MenuTitle,Buttons);
    else
      k=menu(MenuTitle,Buttons);
    end

    if ~isempty(k),
      if k >= 6,
        NeighborhoodChanged = 1;
      end

      switch k,
      case 1,
        Display(1).colorstate = 1 - Display(1).colorstate;
      case 2,
        Display(1).displaybackbone = 1 - Display(1).displaybackbone;
      case 3,
        Display(1).labelbases = fontsize - Display(1).labelbases;
      case 4,
        Display(1).superimposestate = Display(1).superimposestate + 1;
        if Display(1).superimposestate > 2,
          Display(1).superimposestate = 0;
        end
        switch Display(1).superimposestate,
        case 0,
          Display(1).superimpose = 0;
          Display(1).superbackbone = 0;
        case 1,
          Display(1).superimpose = 1;
          Display(1).superbackbone = 1;
        case 2,
          Display(1).superimpose = 1;
          Display(1).superbackbone = 0;
        end
      case 5,
        Display(1).nearbyatoms = 1 - Display(1).nearbyatoms;
      case 6,
        Display(1).neighborhood = zeros(1,8);
      case 7,
        Display(1).neighborhood(1) = 1;
      case 8,
        Display(1).neighborhood(2) = 1;
        Display(1).neighborhood(3) = 0;
        Display(1).neighborhood(4) = 0;
      case 9,
        Display(1).neighborhood(2) = 0;
        Display(1).neighborhood(3) = 1;
        Display(1).neighborhood(4) = 0;
      case 10,
        Display(1).neighborhood(2) = 0;
        Display(1).neighborhood(3) = 0;
        Display(1).neighborhood(4) = 1;
      case 11,
        Display(1).neighborhood(5) = 1;
      case 12,
        Display(1).neighborhood(6) = 1;
      case 13,
        Display(1).neighborhood(7) = 1;
      case 14,
        Display(1).neighborhood(8) = 8;
      case 15,
        Display(1).neighborhood(8) = 12;
      case 16,
        Display(1).neighborhood(8) = 16;
      case 17,
        Display(1).backbonetrace = 1 - Display(1).backbonetrace;
        Display(1).aabackbonetrace = 1 - Display(1).aabackbonetrace;
      end
    end
% end
end
