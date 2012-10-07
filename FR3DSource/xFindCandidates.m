% xFindCandidates finds and ranks candidate fits to Model

function [Found,PS] = xFindCandidates(File,Model,Verbose)

Found = uint16([]);                         % store indices as integers

starttime = cputime;                        % to track time usage

NNZ = [];                                   % number of entries in PS

for f=1:length(File),
 if File(f).NumNT > 0,                      % if the file is non-empty
  filestarttime = cputime;                  % for this file only
  Codes = cat(1,File(f).NT(:).Code);        % codes to use in nucleotide mask

  % the following check is silly, because these fields always exist,
  % but might be empty, which is not being checked for

  N = length(Codes);

  %-----------------------------------------------------------------------
  % go through all pairs and set up sparse screening matrices
  % This could be sped up by combining the later preprocessing step with
  % this one so that fewer calculations are done after one pairwise screen
  % matrix is found that has few non-zero entries.
  % In that case, the order in which the pairwise screen matrices are set
  % up will matter.  Best to start with a pair that is restricted to a
  % certain basepair or stack or paircode.

  for i=2:Model.NumNT                % loop through model nucleotides
    for j=1:(i-1)
      PS{i,j} = xPairwiseScreen(File(f),Codes,Model,i,j);
      PS{j,i} = PS{i,j}';
      NNZ(i,j) = nnz(PS{i,j});              % number of non-zero entries
      NNZ(j,i) = NNZ(i,j);
    end
  end

  %-----------------------------------------------------------------------

  if 1 > 0,                                 % this preprocessing works

    % this code removes non-zero entries of PS{i,j}
    % by simply comparing PS{i,1}, PS{i,2}, etc. and only keeping
    % those rows which have some non-zero entry somewhere in them.
    % If PS{i,j}(k,:) is all zeros, then k cannot correspond to i,
    % so zero out PS{i,m}(k,:) for all other m.

    % First, order the nucleotides so that the tightest constraint
    % is used first to impact the other constraints:

    Perm = xOrderQueryNucleotides(NNZ);         % find a good ordering
    PS   = PS(Perm,Perm);                       % re-order nucleotides

    NewNNZ = [];

    % Next, determine which rows of PS{i,j} are empty for at least one j

    for i = 2:Model.NumNT,
      RowSums = zeros(N,1);
      for j = 1:Model.NumNT,
        if j ~= i,
          RowSums = RowSums + any(PS{i,j},2);
        end
      end

      % Identify the "good" and "bad" rows; both are useful

      badi  = find(RowSums < Model.NumNT-1);
      goodi = find(RowSums == Model.NumNT-1);

      % Two ways to zero out, depending on whether just a few rows
      % are being zeroed out, or many.

      if length(badi) < N/2,
        for j = 1:Model.NumNT,
          if j ~= i,
            PS{i,j}(badi,:) = 0 * PS{i,j}(badi,:);
            PS{j,i} = PS{i,j}';
          end
        end
      elseif length(goodi) <= N/2,
        for j = 1:Model.NumNT,
          if j ~= i,
            A = sparse(N,N);
            A(goodi,:) = PS{i,j}(goodi,:);
            PS{i,j} = A;
            PS{j,i} = A';
          end
        end
      end
    end

    PS(Perm,Perm) = PS;             % undo the permutations

    % Recalculate the number of non-zero entries

    for i = 1:Model.NumNT,
      for j = (i+1):Model.NumNT,
        a = nnz(PS{i,j});
        NewNNZ(i,j) = a;
        NewNNZ(j,i) = a;
      end
    end

    NNZ = NewNNZ;
  end

  %-----------------------------------------------------------------------

  Perm = xOrderQueryNucleotides(NNZ);         % find a good ordering
  PS   = PS(Perm,Perm);                       % re-order nucleotides

  NN = min(9,Model.NumNT);                    % start with at most 9

  [List,SS] = xFindPolyhedra(Model,NN,PS(1:NN,1:NN));

  %-----------------------------------------------------------------------

  for r = 10:Model.NumNT,              % additional nucleotides, if needed
    if Verbose > 1,
      tic
    end
    if length(List(:,1)) > 0,
      [List,SS] = xAddNucleotide(Model,List,PS,SS,r);
    end                                % end if length(List(:,1)) > 0,
    if Verbose > 1,
      fprintf('Adding nucleotide %2d took    %14.3f seconds\n',r,toc);
    end
  end                                  % end for r

  %-------------------------------------------------------------------------

  if ~isempty(List),
    List(:,Perm) = List(:,1:Model.NumNT);     % re-order nucleotides
  end

  [s,t] = size(List);
  if s > 0,
    Found = [Found; [List uint16(f*ones(s,1))]];
    if Verbose > 0,
      fprintf('Found %7d possibilities from %10s in %8.3f seconds\n', s, File(f).Filename, cputime-filestarttime);
    end
  end

  drawnow

 end % if File(f).NumNT > 0
end  % end of for loop length(File)

[s,t] = size(Found);

if (length(File) > 1) && (Verbose > 0),
  fprintf('Found %7d possible candidates in %8.3f seconds\n', s, (cputime-starttime));
end

drawnow

PS(Perm,Perm) = PS;             % undo the permutations
