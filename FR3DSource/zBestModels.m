% zBestModels(File) determines the model or models File that has the "best" 3D structure according to some criteria.  If Verbose > 0, it writes to the current file fid.  This will be useful for NMR files.

function [Models] = zBestModels(File,Verbose)

if nargin < 2,
  Verbose = 0;
end

for f = 1:length(File);
  if length(File(f).NT) > 0,
    clear ModelInd ModelRed
    Model = cat(2,File(f).NT.ModelNum);                  % models
    U = unique(Model);
    for u = 1:length(U),                              % loop through chains
      ModelInd{u} = find(Model == U(u));              % NTs in this chain
    end
                             % tally redundancies between all pairs of chains
    for u = 1:length(U),
      for v = 1:length(U),
        ModelRed(u,v) = sum(sum(File(f).Redundant(ModelInd{u},ModelInd{v})));
      end
    end

    KeepModels = [];

    if sum(sum(ModelRed)) > 0,                        % redundant chains exist
      clear BPQual IntQual
      E = File(f).Edge;
      Q = (abs(E) > 0) .* (abs(E) < 13);              % indicator of edges
      R = (abs(E) > 0) .* (abs(E) < 25);              % indicator of edges
      for u = 1:length(U),
        for v = 1:length(U),
          BPQual(u,v) = sum(sum(Q(ModelInd{u},ModelInd{v})));
          IntQual(u,v) = sum(sum(R(ModelInd{u},ModelInd{v})));
        end
        BPQual(u,u) = 1;                              % in case no BP w/in chain
      end

      % extend by transitivity to find chains that are connected by basepairs

      QualP = full((BPQual + BPQual^2 + BPQual^3 + BPQual^4 + BPQual^5) > 0);

      if Verbose > 1,
        clf
        U
        zCircularDiagram(File(f))
        drawnow
        full(ModelRed)
        full(IntQual)
        QualP
      end

      if all(QualP > 0),
        KeepModels = 1:length(U);
        if Verbose > 1,
          fprintf('Keeping all models\n');
        end
      else
        G  = 1-(sum(QualP) > 0);             % chains used, or no inter
        m  = 0;                              % maximum quality so far
        n  = 0;
        while any(G == 0),                   % unused chains
          i = min(find(G == 0));              % first unused chain
          KC = find(QualP(i,:) > 0);          % chains with 1st chain
          G(KC) = ones(1,length(KC));          % these have been considered
          Q  = sum(sum(BPQual(KC,KC)));       % quality of these chains
          QInt = sum(sum(IntQual(KC,KC)));
          if Q > m,
            KeepModels = KC;
            m = Q;
            n = QInt;
          elseif Q == m,
            if QInt > n,
              KeepModels = KC;
              n = QInt;
            end
          end
        end
        if Verbose > 1,
          fprintf('Keeping model(s) ');
          fprintf('%d ', KeepModels);
          fprintf('%s ', U(KeepModels));
          fprintf('\n');
        end
      end

      if Verbose > 1,
        pause
      end
    else
      KeepModels = 1:length(U);                 % keep them all
    end

    Models{f} = U(KeepModels);

    if Verbose > 0,
      fprintf('%s\t',File(f).Filename);
      fprintf('%s',U(KeepModels));
      fprintf('\n');
    end

  else
    Models{f} = '';
  end
end
