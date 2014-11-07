% pConsensusInteractions(Search) determines the consensus interactions in the candidates in Search

function [Edge,BPh,BR,Search] = pConsensusInteractions(Search,Verbose)

if nargin < 2,
  Verbose = 0;
end

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

f = Search.Candidates(:,N+1);           % file numbers of motifs

% ----- Calculate coplanar measure for each File in Search

if ~isfield(Search.File(1),'Coplanar'),

  disp('Adding coplanar values')

  clear NewFile
  for ff = 1:length(Search.File),
    F = Search.File(ff);
    if ~isempty(F.NT),
      F.Coplanar = sparse(F.NumNT,F.NumNT);
      NewFile(ff) = F;
    end
  end
  Search.File = NewFile;

  for c = 1:length(Search.Candidates(:,1)),
    ff = Search.Candidates(c,N+1);
    i = Search.Candidates(c,1:N);
    for a = 1:length(i),
      for b = (a+1):length(i),
        if Search.File(ff).Edge(i(a),i(b)) ~= 0,
          NT1 = Search.File(ff).NT(i(a));
          NT2 = Search.File(ff).NT(i(b));
          Pair = zAnalyzePair(NT1,NT2);
          Search.File(ff).Coplanar(i(a),i(b)) = Pair.Coplanar;
          Search.File(ff).Coplanar(i(b),i(a)) = Pair.Coplanar;
        end
      end 
    end
  end
end

% -------------------------------------------- find consensus pairs and stacks

Edge = sparse(N,N);                             % place to store consensus

for a = 1:N,                                    % first NT of possible pair
  for b = (a+1):N,                              % second NT of possible pair
    e = [];                                     % record observed edges
    cp= [];                                     % coplanar, for basepairs
    for c = 1:L,                                % run through candidates
      i = Search.Candidates(c,a);               % index of first nucleotide
      j = Search.Candidates(c,b);               % index of second nucleotide
      e = [e Search.File(f(c)).Edge(i,j)];      % append observed interaction
      
      % modified by Anton to deal with old motif atlas releases 12/6/2011
      if isfield(Search.File(f(c)),'LooseCoplanar')
          cp= [cp Search.File(f(c)).LooseCoplanar(i,j)]; % append coplanar information
      else
          cp= [cp 0];
      end
     % modified by Anton to deal with old motif atlas releases 12/6/2011

    end

    e = fix(e);                                 % round subcategories

    for d = 1:length(e),                        % loop through instances
      if any(e(d) == [-1 -2 -7 -8 -12 -22 -23 -101 -102 -107 -108 -112 -122 -123]),
        e(d) = -e(d);                           % don't distinguish sign here
                                                % these are symmetric enough
      end
    end

    ae = (abs(e)>0) .* (abs(e) < 30) .* e;      % FR3D basepairs and stacks

    anc= (abs(e)>100) .* (abs(e) < 120) .* (cp > 0) .* e; % npair but coplanar

    % ------------------------------------------- determine the "best" inter

    if length(ae(ae~=0)) > 0,                   % if there is a true interact
      emode = mode(ae(ae~=0));                  % store the code for it
    elseif length(anc(anc~=0)) > 0,             % if there is a near coplanar
      emode = mode(anc(anc~=0));                % store the code for it
      emode = sign(emode)*(abs(emode)-100);     % corresponding true code
    else
      emode = NaN;
    end

    nemode = sign(emode) * (abs(emode) + 100);  % near code of mode

    numeqemode = sum(emode == e);               % number equal to the mode

    numeqnemode = sum(nemode == e);             % number equal to near mode

    numeqnemodec= sum(nemode == anc);           % coplanar and equal near mode

    T = numeqemode;                             % true
    NC = numeqnemodec;                          % near but coplanar
    NN = numeqnemode - NC;                      % just plain near

    if Verbose > 0,
      figure(57)
      xlabel('True and near coplanar pairs over length');
      ylabel('Near basepairs over length')
      hold on
      axis([-0.1 1.1 -0.1 1.1]);
    end

    if 3*T + 2*NC + NN > L && (T+NC) > L/8 || 3*T + 2*NC > 30,       % new criterion
      Edge(a,b) = emode;
      if Verbose > 0,
        if Edge(a,b) == 0 && abs(emode) < 20,
          fprintf('pConsensusInteractions: +++++++++++++++++++++ New basepair consensus\n');
          fprintf('pConsensusInteractions: Bases %d and %d\n', a, b);
          fprintf('pConsensusInteractions: Consensus used to be %s\n', zEdgeText(Edge(a,b)));
          fprintf('pConsensusInteractions: Consensus is now %s\n', zEdgeText(emode));
          fprintf('pConsensusInteractions: %d true, %d near coplanar, and %d near out of %d instances\n',full(T),full(numeqnemodec),full(NN),L);
          plot((T+NC)/L+0.01*randn,NN/L+0.01*randn,'g.');
        end
        if abs(Edge(a,b)) < 20,
          plot((T+NC)/L+0.01*randn,NN/L+0.01*randn,'k.');
        end
      end
    else
      Edge(a,b) = 0;
      if Verbose > 0,
        if Edge(a,b) ~= 0 && abs(Edge(a,b)<20),
          fprintf('pConsensusInteractions: --------------------- Lost basepair consensus\n');
          fprintf('pConsensusInteractions: Bases %d and %d\n', a, b);
          fprintf('pConsensusInteractions: Consensus used to be %s\n', zEdgeText(Edge(a,b)));
          fprintf('pConsensusInteractions: %d true, %d near coplanar, and %d near out of %d instances\n',full(T),full(numeqnemodec),full(NN),L);
          plot((T+NC)/L+0.01*randn,NN/L+0.01*randn,'r.');
  %        pause
        end
      end
    end

    Edge(b,a) = -Edge(a,b);                     % record twice

  end
end

% ----------------------------------------------- BPh and BR edges used

% BBEdge{1} = 

% ----------------------------------------------- consensus BPh interactions

% base-backbone edges for illustration
edgemap = {'B',' W',' H',' H',' W',' H',' H',' H',' H',' S',' W',' W',' W',' H',' W',' H',' H',' H',' W'};
edgemap = {[1 2], 2, 3, 0, 2, 3, 3, 3, 0, 1, 2, 2, 2, 0, 2, 3, 0, 3, 2};
edgename = 'SWH';
% 0 is 0BPh/BR, 1 is sugar, 2 is WC, 3 is Hoogsteen
% Note:  0BPh and 0BR are stored as 0, which is the same as no interaction

BPh = zeros(N,N);

for a = 1:N,                                      % first NT of possible BPh
  for b = 1:N,                                    % second NT of possible BPh
%    if a ~= b,                                    % don't bother with self interactions
      edge = [];                                  % record observed edges
      nearedge = [];                              % record edge of near interaction
      for c = 1:L,                                % run through candidates
        i = Search.Candidates(c,a);               % index of first nucleotide
        j = Search.Candidates(c,b);               % index of second nucleotide
        e = Search.File(f(c)).BasePhosphate(i,j); % append observed interaction

        if (e > 0) && (e < 20),
          edge = [edge edgemap{e}];               % append edge of interaction
        elseif (e > 100 && e < 120),              % near interaction
          nearedge = [nearedge edgemap{e-100}];   % append edge of near interaction
        end
      end

      if length(edge) > 0,                        % at least one interaction
        [truemode,numtruemode] = mode(edge);
        numneartruemode        = sum(nearedge == truemode);

        if numtruemode >= L/4 && 2*numtruemode + numneartruemode > L && truemode > 0,
          BPh(a,b) = truemode;                    % use this edge
        end

        if Verbose > 0 && numtruemode > 0 && a == b && truemode > 0,
          fprintf('pConsensusInteractions: BPh interaction between %d and %d, %d candidates\n', a, b, L);
          edge
          nearedge
          if BPh(a,b) > 0,
            fprintf('pConsensusInteractions: BPh interaction on %s edge between %d and %d, %d candidates\n', edgename(truemode), a, b, L);
          else
            fprintf('pConsensusInteractions: No agreement on BPh\n');
          end

%          pause
        end
      end
%    end
  end
end

% ----------------------------------------------- consensus BR interactions

if ~isfield(Search.File(1),'BaseRibose'),
  disp('Adding base-ribose interactions')
  Search.File = zBaseRiboseInteractions(Search.File);
end

BR  = zeros(N,N);

for a = 1:N,                                    % first NT of possible BPh
  for b = 1:N,                                  % second NT of possible BPh
%    if a ~= b,                                  % don't bother with self interactions
      edge = [];                                  % record observed edges
      nearedge = [];                              % record near interaction
      for c = 1:L,                                % run through candidates
        i = Search.Candidates(c,a);               % index of first nucleotide
        j = Search.Candidates(c,b);               % index of second nucleotide
        e = Search.File(f(c)).BaseRibose(i,j);    % append observed interaction

        if (e > 0) && (e < 20),
          edge = [edge edgemap{e}];               % append actual interaction
        elseif (e > 100 && e < 120),              % near interaction
          nearedge = [nearedge edgemap{e-100}];   % append near interaction
        end
      end

      if length(edge) > 0,                        % at least one interaction
        [truemode,numtruemode] = mode(edge);
        numneartruemode        = sum(nearedge == truemode);

        if numtruemode >= L/4 && 2*numtruemode + numneartruemode > L && truemode > 0,
          BR(a,b) = truemode;                      % use this edge
        end

        if Verbose > 0 && numtruemode > 0 && a == b && truemode > 0,
          fprintf('pConsensusInteractions: BR interaction between %d and %d, %d candidates\n', a, b, L);
          edge
          nearedge
          if BR(a,b) > 0,
            fprintf('pConsensusInteractions: BR interaction on %s edge between %d and %d, %d candidates\n', edgename(truemode), a, b, L);
          else
            fprintf('pConsensusInteractions: No agreement on BR\n');
          end

%          pause
        end
      end
%    end
  end
end

if 0 > 1,
  % ----------------------------------------------- consensus BB interactions?
  % ------------ None found in IL/1.6

  BB  = zeros(N,N);

  for a = 1:N,                                    % first NT of possible BPh
    for b = 1:N,                                  % second NT of possible BPh
%      if a ~= b,                                  % don't bother with self interactions
        edge = [];                                  % record observed edges
        nearedge = [];                              % record near interaction
        for c = 1:L,                                % run through candidates
          i = Search.Candidates(c,a);               % index of first nucleotide
          j = Search.Candidates(c,b);               % index of second nucleotide
          e = Search.File(f(c)).BasePhosphate(i,j);    % append observed interaction

          if e == 0,
            e = Search.File(f(c)).BaseRibose(i,j);  % use BR too
          end

          if e == 0 && j > 1,
            e = Search.File(f(c)).BaseRibose(i,j-1); % use BR from neighboring base too
          end

          if (e > 0) && (e < 20),
            edge = [edge edgemap{e}];               % append actual interaction
          elseif (e > 100 && e < 120),              % near interaction
            nearedge = [nearedge edgemap{e-100}];   % append near interaction
          end

        end

        if length(edge) > 0,                        % at least one interaction
          [truemode,numtruemode] = mode(edge);
          numneartruemode        = sum(nearedge == truemode);

          if numtruemode >= L/4 && 2*numtruemode + numneartruemode > L && truemode > 0,
            BR(a,b) = truemode;                      % use this edge
          end

          if BB(a,b) > 0 && BPh(a,b) == 0 && BR(a,b) == 0,  % BB but not BPh or BR --> so something new!

            fprintf('pConsensusInteractions: BB interaction between %d and %d, %d candidates\n', a, b, L);
            edge
            nearedge
            if BB(a,b) > 0,
              fprintf('pConsensusInteractions: BB interaction on %s edge between %d and %d, %d candidates\n', edgename(truemode), a, b, L);
            end

%            pause
          end
        end
      end
%    end
  end
end
