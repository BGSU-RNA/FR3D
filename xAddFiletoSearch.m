% xAddFiletoSearch(File,Search) adds data on nucleotides found in Search
% It assumes that File is in the proper order for the index in Candidates

function [Search] = xAddFiletoSearch(File,Search)

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

if (s==0),
  Search.CandidateFilenames{1} = '';
  Search.File(1).Filename = '';

else

for i = 1:s,
  f = double(Candidates(i,N+1));        % file number for this candidate
  Search.CandidateFilenames{f} = Search.Filenames{f};
end

if ~isempty(File),
  for f = 1:max(Candidates(:,N+1)),
    Search.File(f).Edge = sparse(zeros(1,1));
  end

  for i = 1:s,
    f = double(Candidates(i,N+1));             % file number for this candidate

    Search.CandidateFilenames{f} = File(f).Filename;
    Search.File(f).Filename = File(f).Filename;
    Search.File(f).NumNT    = File(f).NumNT;
    Search.File(f).Info     = File(f).Info;

    Indices = Candidates(i,1:N);               % indices of nucleotides

    for j = Indices,
      Search.File(f).NT(j) = File(f).NT(j);
      for k = Indices,
        Search.File(f).Edge(j,k) = File(f).Edge(j,k);
        Search.File(f).Edge(k,j) = File(f).Edge(k,j);
      end
    end
  end
end

end