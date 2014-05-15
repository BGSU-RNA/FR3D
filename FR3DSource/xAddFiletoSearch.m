% xAddFiletoSearch(File,Search) adds data on nucleotides found in Search.  It assumes that File is in the proper order for the index in Candidates

function [Search] = xAddFiletoSearch(File,Search)

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

if (s==0),                                % no candidates
  Search.CandidateFilenames{1} = '';
  Search.File(1).Filename = '';
else

  [y,p] = sort(double(Candidates(1,1:N))); % put nucleotides of first candidate in increasing order
  
  Cand = double(Candidates(:,p));     % re-order nucleotides acc to first

  if isfield(Query,'MaxDiffMat'),
    MaxDiff = diag(Query.MaxDiffMat(p,p),1);
  else
    MaxDiff = Inf*ones(1,N-1);
  end
  
  % ---------------------------- Calculate maximum gaps between cand. nucleotides
  
  maxinsert = zeros(1,N-1);
  for c = 1:s,
    maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
  end
  
  % ---------------------------- Add nucleotide information
  
  if ~isempty(File),
    for f = 1:max(Candidates(:,N+1)),            % loop through files
      Search.File(f).Edge     = sparse([],[],[],N,N);
      Search.File(f).BasePhosphate = sparse([],[],[],N,N);
      Search.File(f).BaseRibose = sparse([],[],[],N,N);
      Search.File(f).Range    = sparse([],[],[],N,N);
      Search.File(f).Crossing = sparse([],[],[],N,N);
      Search.File(f).Covalent = sparse([],[],[],N,N);
      Search.File(f).Backbone = sparse([],[],[],N,N);
      Search.File(f).Coplanar = sparse([],[],[],N,N);
      Search.CandidateFilenames{f} = File(f).Filename;
      Search.File(f).Filename      = File(f).Filename;
      Search.File(f).PDBFilename   = File(f).PDBFilename;
      Search.File(f).NumNT         = File(f).NumNT; % max number, some empty
      Search.File(f).Info          = File(f).Info;

      AllIndices{f} = [];
      AllOrigIndices{f} = [];
      AllAA{f} = [];
    end
  
  tic

    for i = 1:s,                           % loop through candidates
      f = double(Candidates(i,N+1));       % file number for this candidate
  
      Indices = Cand(i,1:N);               % indices of nucleotides

      Indices = [Indices xNeighborhood(File(f),Indices,[1 0 0 1 1 1 1 16])];  % largest neighborhood

      AllIndices{f} = [AllIndices{f} Indices]; % accumulate indices needed for file f
      AllOrigIndices{f} = [AllOrigIndices{f} Cand(i,1:N)];

      % include intervening nucleotides, if only a few, for alignments
  
      for n = 1:(N-1),                               % loop through nucleotides
        if (MaxDiff(n) < Inf) || (maxinsert(n) < 5),   % if only few insertions
          if double(Indices(n+1)) - double(Indices(n)) > 1,   % increasing order
            AllIndices{f} = [AllIndices{f} (Indices(n)+1):(Indices(n+1)-1)];
          elseif double(Indices(n+1)) - double(Indices(n)) < -1, % dec order
            AllIndices{f} = [AllIndices{f} (Indices(n)-1):-1:(Indices(n+1)+1)];
          end
        end
      end
    end

fprintf('Finding nearby nucleotides, ');
toc
zFlushOutput
tic

    for f = 1:length(AllIndices),
      j = unique(AllIndices{f});        % only do each index once
      Search.File(f).NT(j) = File(f).NT(j);
      Search.File(f).Edge(j,j) = File(f).Edge(j,j);
      Search.File(f).Edge(j,j) = File(f).Edge(j,j);
      Search.File(f).BasePhosphate(j,j) = File(f).BasePhosphate(j,j);
      Search.File(f).BasePhosphate(j,j) = File(f).BasePhosphate(j,j);
      Search.File(f).BaseRibose(j,j) = File(f).BaseRibose(j,j);
      Search.File(f).BaseRibose(j,j) = File(f).BaseRibose(j,j);
      Search.File(f).Range(j,j) = File(f).Range(j,j);
      Search.File(f).Range(j,j) = File(f).Range(j,j);
      Search.File(f).Covalent(j,j) = File(f).Covalent(j,j);
      Search.File(f).Covalent(j,j) = File(f).Covalent(j,j);
      Search.File(f).Crossing(j,j) = File(f).Crossing(j,j);
      Search.File(f).Crossing(j,j) = File(f).Crossing(j,j);
      Search.File(f).Backbone(j,j) = File(f).Backbone(j,j);
      Search.File(f).Backbone(j,j) = File(f).Backbone(j,j);
      Search.File(f).Coplanar(j,j) = File(f).Coplanar(j,j);
      Search.File(f).Coplanar(j,j) = File(f).Coplanar(j,j);
    end

fprintf('Adding nucleotides and interaction information, ');
toc
zFlushOutput

    for f = 1:length(AllOrigIndices),
      AllOrigIndices{f} = unique(AllOrigIndices{f});
      if isfield(File(f),'AA') && length(File(f).AA) > 0,
        c = cat(1,File(f).NT(AllOrigIndices{f}).Center);      % nucleotide centers
        a = cat(1,File(f).AA.Center);               % amino acid centers
        [ss,tt] = size(a);
        [sss,ttt] = size(c);
        if ss == 0 || sss == 0,
          Search.File(f).AA = [];
        else
          D = zDistance(c,a);                         % NT-AA distances
          [i,j,k] = find(D);                          % 
          w = find(k < 16);                           % within 16 Angstroms
          u = unique(j(w));                           % indices of amino acids
          Search.File(f).AA(u) = File(f).AA(u);       % add these AAs
        end
      else
        Search.File(f).AA = [];
      end
    end

fprintf('Finding nearby amino acids, ');
toc
zFlushOutput

  elseif ~isfield(Search,'CandidateFilenames'),

    List = {};
    fprintf('Attempting to convert filenames and lists to current format\n');
    for j=1:length(Search.Filenames),
      if strfind(Search.Filenames{j},'_list'),
        fprintf('If %s has changed or is no longer present, this may fail!\n',Search.Filenames{j});
      end
      List = [List; zReadPDBList(Search.Filenames{j},1)];
    end
    Search.CandidateFilenames = List;
  end
end
