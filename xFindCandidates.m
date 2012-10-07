% xFindCandidates finds and ranks candidate fits to Model

function Found = xFindCandidates(File,Model)

Found = uint16([]);                         % store indices as integers

starttime = cputime;                        % to track time usage

for f=1:length(File),
 if File(f).NumNT > 0,                      % if the file is non-empty
  filestarttime = cputime;                  % for this file only
  Codes = cat(1,File(f).NT(:).Code);        % codes to use in nucleotide mask
  if isfield(Model,'OKPairs') | isfield(Model,'ExPairs'),%  screen by paircode
    N = length(Codes);
    PC = Codes * ones(1,N) + 4*(ones(N,1)*(Codes'-1));
  else
    PC = [];
  end

  for i=2:min(9,Model.NumNT)                % loop through model nucleotides
    for j=1:(i-1)
      PS{i,j} = xPairwiseScreen(File(f),Codes,Model,i,j,PC);
      PS{j,i} = PS{i,j}';
      NNZ(i,j) = nnz(PS{i,j});              % number of non-zero entries
      NNZ(j,i) = NNZ(i,j);
    end
  end

  Perm = sFindPermutation(NNZ);               % find best nucleotide order
  PS   = PS(Perm,Perm);                       % re-order nucleotides

  [List,SS] = xFindPolyhedra(Model,min(9,Model.NumNT),PS);

  if ~isempty(List),
    List(:,Perm) = List(:,1:min(9,Model.NumNT));     % re-order nucleotides
  end

  %-----------------------------------------------------------------------
  for r = 10:Model.NumNT,              % additional nucleotides, if needed
    if length(List(:,1)) > 0,
      for q = 1:(r-1),
        S{r,q} = xPairwiseScreen(File(f),Codes,Model,r,q);
      end                              % end for q
        [List,SS] = xAddNucleotide(Model,List,S,SS,r);
    end                                % end if length(List(:,1)) > 0,
%   fprintf('Adding nucleotide %2d took     %14.3f seconds\n',r,toc);
  end                                  % end for r

  %-------------------------------------------------------------------------

  [s,t] = size(List);
  if s > 0,
    Found = [Found; [List uint16(f*ones(s,1))]];
    fprintf('Found %7d possibilities from %10s in %8.3f seconds\n', s, File(f).Filename, cputime-filestarttime);
  end

  drawnow

 end % if File(f).NumNT > 0
end  % end of for loop length(File)

[s,t] = size(Found);
if length(File) > 1,
  fprintf('Found %7d possible candidates in %8.3f seconds\n', s, (cputime-starttime));
end

drawnow

