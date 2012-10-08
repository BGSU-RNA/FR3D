% zMotifSignatures(Motif,Strands,Rotation,Type) determines a motif signature for the motif whose FR3D query information is encoded in Motif.  Strands is the number of strands in the motif, Rotation tells which strand to start with, Type is 0 for basepairs only, 1 for basepairs and near pairs, 2 for pairs and stacks

function [Sig,AllSig,AllPhoneme] = zMotifSignature(Edge,Strands,Rotation,Type)

if nargin < 4,
  Type = 0;                                       % basepairs only
end

if nargin < 3,
  Rotation = 1;
end

if nargin < 2,
  Strands = 2;
end

N = length(Edge(1,:));                            % number of nucleotides

E = abs(fix(Edge));                               % abs values of inter codes

if Type == 0,
  E = E .* (E < 14);                              % basepairs only
elseif Type == 1,
  E = E .* (E < 14) + E .* (E > 100) .* (E < 114); % basepairs and near pairs
elseif Type == 2,
  E = E .* (E < 24);                              % basepairs and stacks
end

Sig = '';                                         % signature starts blank
Pho = '';                                         % phoneme starts blank

% ----------------------------------------------- % determine strand boundaries

if Strands == 1,
  S = N-1;
elseif Strands == 2,
  if E(1,N) == 1,                            % flanking cWW pair
    i = N;
    CWW = 0;
    while CWW == 0 && i > 0,
      for j = (i+1):N,
        if E(i,j) == 1,
          CWW = 1;
        end
      end
      if CWW == 0,
        i = i - 1;
      end
    end
                                             % probe for more in this strand
%    while Motif.DifferenceSign(i,i+1) == -1 && i < N-1,
%      i = i + 1;
%    end
    S = i;
  else
    S = N;                                   % go through all rows
  end

  if Rotation == 2;
    p = [(S+1):N 1:S];                       % re-order strands
    Edge = Edge(p,p);
    E    = E(p,p);
    S = N-S;
  end

end

fprintf('zMotifSignature: Rotation = %d, S = %d\n', Rotation, S);

SS = '';                                      % initialize signature
MRC = N;                                      % maximum remaining column
%full(E)

if Strands == 1,
  for i = 1:S,                                  % go in sequential order
    if sum(E(i,:)) == 0,                        % i makes no interaction
      SS = [SS '-F'];                           % unpaired base on the left
    end

    for j = N:-1:(i+1),                       % work back through columns
      if E(i,j) > 0,
        SS = [SS '-' zEdgeText(Edge(i,j))]; 
      end
    end
  end
elseif Strands == 2,
  for i = 1:S,
    if sum(E(i,:)) == 0,                        % i makes no interaction
      SS = [SS '-L'];                           % unpaired base on the left
    end

    while sum(E(i:S,MRC)) == 0,                 % no remaining interactions in MRC
      if sum(E(:,MRC)) == 0,                    % NT in MRC makes no interactions
        SS = [SS '-R'];
      else                                      % look for interaction on right
        for k = (S+1):(MRC-1),                  % loop through rows
          if E(k,MRC) > 0,
            SS = [SS '-' zEdgeText(Edge(k,MRC))]; 
          end
        end
      end
      MRC = MRC - 1;                            % move left
    end

    for j = N:-1:(i+1),                       % work back through columns
      if E(i,j) > 0,
        SS = [SS '-' zEdgeText(Edge(i,j))]; 
      end

      if i == S && MRC > j && j > S+1,
        MRC = j;                              % need to check this column
        if sum(E(:,MRC)) == 0,                % NT in MRC makes no interactions
          SS = [SS '-R'];
        else                                  % look for interaction on right
          for k = (S+1):(MRC-1),              % loop through rows
            if E(k,MRC) > 0,
              SS = [SS '-' zEdgeText(Edge(k,MRC))]; 
            end
          end
        end
      end
    end
  end
end

SS = strrep(SS,' ','');
while SS(1) == '-' && length(SS) > 1,
  SS = SS(2:end);
end

while SS(end) == '-' && length(SS) > 1,
  SS = SS(1:(end-1));
end

Sig = SS;

fprintf('zMotifSignature: Signature: %s\n', SS);

% ----------------------------------------- clean up signature

while Sig(1) == '-' && length(Sig) > 1,
  Sig = Sig(2:end);
end

while Sig(end) == '-' && length(Sig) > 1,
  Sig = Sig(1:(end-1));
end

Pho = Sig;
Phoneme = zPhoneme;
for i = 1:length(Phoneme(:,1)),
  Pho = strrep(Pho,Phoneme{i,1},Phoneme{i,2});
end
Pho = strrep(Pho,'-','');
Pho = lower(Pho);

fprintf('zMotifSignature: Phonetic:  %s\n', Pho);

AllSig{Rotation} = Sig;
AllPhoneme{Rotation} = Pho;

if Strands > 1 && Rotation < Strands,
  [Si,Al,Ph] = zMotifSignature(Edge,Strands,Rotation+1,Type);
  AllSig{Rotation+1} = Si;
  AllPhoneme{Rotation+1} = Ph{Rotation+1};
end

