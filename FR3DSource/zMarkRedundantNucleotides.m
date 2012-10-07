% zMarkRedundantNucleotides(File) identifies redundant chains in File and sets File.Redundant(i,j) = 1 if nucleotides with indices i and j, from different chains, are redundant when these chains are aligned and superimposed.

% Some chains are identified with different letters.
% In NMR structures, however, they may all be called the same thing

function [File] = zMarkRedundantNucleotides(File, Verbose)

if nargin < 2,
  Verbose = 1;
end

DiscCutoff = 0.4;                               % limit on how dissimilar

for f = 1:length(File),

 File(f).Redundant = sparse([],[],[],length(File(f).NT),length(File(f).NT));
 LongestChain{f} = [1 1];

 ChainIndicator = zeros(1,length(File(f).NT));
 clear NuclNum

 if length(File(f).NT) > 0,                     % File has nucleotides
  Chain = cat(2,File(f).NT.Chain);              % all chain identifiers
  U = unique(Chain);                            % unique chain identifiers

  for u = 1:length(U),                          % loop through all chains IDs
    i = find(Chain == U(u));                    % indices of this chain ID
    ChainIndicator(i) = u;                      % mark these as all the same
    
    for b = 1:length(i),
      NuclNum{b} = File(f).NT(i(b)).Number;     % store nucleotide numbers
    end

    if length(unique(NuclNum)) < length(i)/2,   % numbers are being re-used
      Connected = sparse(eye(length(i)));
      for b = 1:(length(i)-1),
        if File(f).Covalent(i(b),i(b+1)) > 0,
          Connected(b,b+1) = 1;                 % covalent connection
          Connected(b+1,b) = 1;
        end
      end
      for b = 1:log2(length(i)),
        Connected = Connected * Connected;      % extend by transitivity
      end
      if Verbose > 1,
        spy(Connected)
        title([File(f).Filename ' connected components in same chain']);
      end

      for b = 1:length(i),                      % mark distinct chains
        ChainIndicator(i(b)) = u + min(find(Connected(b,:)>0))/10000;
      end
    end
  end

  U = unique(ChainIndicator);                   % unique chain identifiers

  if length(U) > 1,                             % more than one chain
    leng = [];
    clear bases indic

    for u = 1:length(U),                        % loop through chains
      i = find(ChainIndicator == U(u));         % indices of this chain
      bases{u} = cat(2,File(f).NT(i).Base);     % store bases of this chain
      indic{u} = i;                             % store indices of this chain
      leng(u)  = length(i);                     % store length of this chain
    end

    if Verbose > 0,
      fprintf('%s %s\n', File(f).Filename, File(f).Info.ExpTechnique);
      for u = 1:length(U),
        fprintf('%s Chain %s (%5.4f):  %s\n', File(f).Filename, File(f).NT(indic{u}(1)).Chain, U(u), bases{u});
      end
    end

    [y,i] = max(leng);
    LongestChain{f} = [min(indic{i(1)}) max(indic{i(1)})]; 
                                                % indices of longest chain

    % ------------------------------------------- compare chains

    if Verbose > 1,
      fprintf('%s chains with redundancy: ', File(f).Filename);
    end

Discrepancies = zeros(length(U));

    for u = 1:length(U),                        % loop through chains
      for v = (u+1):length(U),                  % loop through chains

        % --------------------------------------- align sequences
        [matches,a,b,ss,tt] = dNeedlemanWunsch(bases{u},bases{v});
        e = find(bases{u}(a) == bases{v}(b));   % locations of agreement
        matches = length(e);                    % number of agreements

        if matches > 1,                         % if there are any matches
          a = a(e);                             % focus on the matches
          b = b(e);

          [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a),File(f),indic{v}(b));

          c = find(abs(A) <= 0.8);              % reasonably similar bases

          if Verbose > 1,
            fprintf('\n%s\n%s\n\n', bases{u}(a), bases{v}(b));
            fprintf('%s\n%s\n', bases{u}(a(c)), bases{v}(b(c)));
          end

          if length(c) > 2 && length(c) > length(a)/2, 
            [Disc,R,MM,CM,A] = xDiscrepancy(File(f),indic{u}(a(c)),File(f),indic{v}(b(c)));

Discrepancies(u,v) = Disc;

            E = eye(length(c),length(c));

            if Disc <= DiscCutoff,
              File(f).Redundant(indic{u}(a(c)),indic{v}(b(c))) = E;
              File(f).Redundant(indic{v}(b(c)),indic{u}(a(c))) = E;

              if Verbose > 1,
                fprintf('%s (%5.4f) %s (%5.4f) ', File(f).NT(indic{u}(1)).Chain, U(u), File(f).NT(indic{v}(1)).Chain, U(v));
              end

            end
          end
        end
      end
    end

    if Verbose > 1,
      fprintf('\n');
    end

    if Verbose > 2,
      figure(2)
      clf
      zCircularDiagram(File(f),1);
    end
    if Verbose > 1,
      figure(1)
      clf
      spy(File(f).Redundant)
      drawnow
    end

  else
    LongestChain{f} = [1 length(File(f).NT)];              % only one chain
  end
 end

 File(f).LongestChain = LongestChain{f};

end

