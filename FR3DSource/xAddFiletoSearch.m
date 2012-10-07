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

    end
  
    for i = 1:s,
      f = double(Candidates(i,N+1));             % file number for this candidate
  
      Search.CandidateFilenames{f} = File(f).Filename;
      Search.File(f).Filename      = File(f).Filename;
      Search.File(f).NumNT         = File(f).NumNT; % max number, some empty
      Search.File(f).Info          = File(f).Info;

      Indices = Cand(i,1:N);               % indices of nucleotides

      Indices = xNeighborhood(File(f),Indices,4);  % largest neighborhood

      for j = Indices,
        Search.File(f).NT(j) = File(f).NT(j);
        for k = Indices,
          Search.File(f).Edge(j,k) = File(f).Edge(j,k);
          Search.File(f).Edge(k,j) = File(f).Edge(k,j);
          Search.File(f).BasePhosphate(j,k) = File(f).BasePhosphate(j,k);
          Search.File(f).BasePhosphate(k,j) = File(f).BasePhosphate(k,j);
          Search.File(f).BaseRibose(j,k) = File(f).BaseRibose(j,k);
          Search.File(f).BaseRibose(k,j) = File(f).BaseRibose(k,j);
          Search.File(f).Range(j,k) = File(f).Range(j,k);
          Search.File(f).Range(k,j) = File(f).Range(k,j);
          Search.File(f).Covalent(j,k) = File(f).Covalent(j,k);
          Search.File(f).Covalent(k,j) = File(f).Covalent(k,j);
          Search.File(f).Crossing(j,k) = File(f).Crossing(j,k);
          Search.File(f).Crossing(k,j) = File(f).Crossing(k,j);
          Search.File(f).Backbone(j,k) = File(f).Backbone(j,k);
          Search.File(f).Backbone(k,j) = File(f).Backbone(k,j);
          Search.File(f).Coplanar(j,k) = File(f).Coplanar(j,k);
          Search.File(f).Coplanar(k,j) = File(f).Coplanar(k,j);
        end
      end

      % include nearby amino acids, for xDisplayCandidates

      if isfield(File(f),'AA'),                

% File(f).Filename

        c = cat(1,File(f).NT(Indices).Center);      % nucleotide centers
        if length(File(f).AA) > 1,
          a = cat(1,File(f).AA.Center);               % amino acid centers
          D = zDistance(c,a);                         % NT-AA distances
          [i,j,k] = find(D);                          % 
          w = find(k < 8);                            % within 8 Angstroms
          u = unique(j(w));                           % indices of amino acids
          Search.File(f).AA(u) = File(f).AA(u);       % add these AAs
        end
      end

      if isfield(File(f),'Het'),
        Search.File(f).Het = [];
      end

      % include intervening nucleotides, if only a few, for alignments
  
      for n = 1:(N-1),                               % loop through nucleotides
       if (MaxDiff(n) < Inf) | (maxinsert(n) < 5),   % if only few insertions
        if double(Indices(n+1)) - double(Indices(n)) > 1,   % increasing order
          for i = (Indices(n)+1):(Indices(n+1)-1),
            Search.File(f).NT(i) = File(f).NT(i);
          end
        elseif double(Indices(n+1)) - double(Indices(n)) < -1, % dec order
          for i = (Indices(n)-1):-1:(Indices(n+1)+1),
            Search.File(f).NT(i) = File(f).NT(i);
          end
        end
       end
      end
    end

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
