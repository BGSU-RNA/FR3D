% zBestChains(File) determines the chain or chains in File that has the "best" 3D structure according to some criteria.  If Verbose > 0, it writes to the current file fid

function [Chains] = zBestChains(File,Verbose)

if nargin < 2,
  Verbose = 0;
end

for f = 1:length(File);
  if length(File(f).NT) > 0,
    clear ChainInd ChainRed
    Chain = cell(length(File(f).NT),1);               % cell array to hold chain identifiers
    [Chain{:}] = deal(File(f).NT.Chain);              % all chain identifiers
    U = unique(Chain);                                % unique chain identifiers
    for u = 1:length(U),                              % loop through chains
      ChainInd{u} = find(ismember(Chain,U{u}));       % indices of NTs in this chain
    end

                             % tally redundancies between all pairs of chains
    ChainRed = [];

    for u = 1:length(U),
      for v = 1:length(U),
        ChainRed(u,v) = sum(sum(File(f).Redundant(ChainInd{u},ChainInd{v})));
      end
    end

    KeepChains = [];

    if sum(sum(ChainRed)) > 0,                        % redundant chains exist
      clear BPQual IntQual
      E = File(f).Edge;                               % basepairs and stacks
      Q = (abs(E) > 0) .* (abs(E) < 13);              % indicator of edges
      R = (abs(E) > 0) .* (abs(E) < 25);              % indicator of edges
      for u = 1:length(U),
        for v = 1:length(U),
          BPQual(u,v) = sum(sum(Q(ChainInd{u},ChainInd{v})));  % number of BP between u and v
          IntQual(u,v) = sum(sum(R(ChainInd{u},ChainInd{v}))); % number of BP and stacks between u and v
        end
        BPQual(u,u) = 1;                              % in case no BP w/in chain
        % BPQual(u,u) = max(1,BPQual(u,u));           % this would be better! 2015-02-02 CLZ
      end

      % extend by transitivity to find chains that are connected by a chain of up to five basepairs

      QualP = full((BPQual + BPQual^2 + BPQual^3 + BPQual^4 + BPQual^5) > 0);

      if Verbose > 1,
        clf
        U
        zCircularDiagram(File(f))
        drawnow
        full(ChainRed)
        full(IntQual)
        QualP
      end

      if all(QualP > 0),
        KeepChains = 1:length(U);
        if Verbose > 1,
          fprintf('Keeping all chains\n');
        end
      else
        G  = 1-(sum(QualP) > 0);             % chains used, or no inter
        m  = 0;                              % maximum quality so far
        n  = 0;
        while any(G == 0),                   % unused chains
          i = min(find(G == 0));             % first unused chain
          KC = find(QualP(i,:) > 0);         % chains connected with 1st chain
          G(KC) = ones(1,length(KC));        % record that these have been considered
          Q  = sum(sum(BPQual(KC,KC)));      % total number of BP between these chains
          QInt = sum(sum(IntQual(KC,KC)));   % total number of BP and stacks between these chains
          if Q > m,                          % better than what has been found so far
            KeepChains = KC;
            m = Q;
            n = QInt;
          elseif Q == m,                     % tie breaker, use total number of interactions
            if QInt > n,
              KeepChains = KC;
              n = QInt;
            end
          end
        end
        if Verbose > 1,
          fprintf('Keeping chains ');
          fprintf('%d ', KeepChains);
          fprintf('%s ', U(KeepChains));
          fprintf('\n');
        end
      end

      if Verbose > 1,
        pause
      end
    else
      KeepChains = 1:length(U);                 % keep them all
    end

    Chains{f} = U(KeepChains);

    if Verbose > 0,
      fprintf('%s\t',File(f).Filename);
      fprintf('%s',U{KeepChains});
      fprintf('\n');
    end

  else
    Chains{f} = '';
  end
end
