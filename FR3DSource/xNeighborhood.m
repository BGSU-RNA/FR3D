% xNeighborhood(File,Indices,v,MaxDiff,MaxInsert) returns indices "near"
% the given indices, in a way determined by the code v

function [NewIndices] = xNeighborhood(File,Indices,v)

N = length(Indices);                                 % number of nucleotides

NewIndices = Indices;                                % always start here

if v > 0,                                          % add basepairs
  for n = 1:N,
    e = abs(File.Edge(Indices(n),:));              % interactions with n
    j = find( (e > 0) .* (e < 14) );               % basepairing with n
    NewIndices = [NewIndices j];
  end
  for n = (N+1):length(NewIndices),                % add what these pair with
    e = abs(File.Edge(NewIndices(n),:));           % interactions with n
    j = find( (e > 0) .* (e < 14) );               % basepairing with n
    NewIndices = [NewIndices j];
  end
end

if v > 1,                                            % add intervening ones
  for n = 1:N,
    NewIndices = [NewIndices (Indices(n)-1):(Indices(n)+1)];  % Anton 1/14/2011 This line causes problems - it adds
    % non-existent indices to the NewIndices structure by going to -1 and +1 positions
  end
end

if v > 2,

  if ~isfield(File,'Distance'),
    File.Distance = [];
  end

  [s,t] = size(File.Distance);

  if isempty(File.Distance) || s < File.NumNT,
    present = zeros(1,File.NumNT);
    c = zeros(File.NumNT,3);
    for i = 1:File.NumNT,
      if i <= length(File.NumNT),
        if ~isempty(File.NT(i).Center),
          c(i,:) = File.NT(i).Center;
          present(i) = 1;
        end
      end
    end
    j = find(present);
    File.Distance = sparse([],[],[],File.NumNT,File.NumNT);
    File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
  end

  d = [1 1 8 10 12];
  a = zeros(1,File.NumNT);
  for j=1:length(Indices),
try
    a = a + (File.Distance(Indices(j),:) < d(v)) .* ...
            (File.Distance(Indices(j),:) > 0);
catch
    fprintf('xNeighborhood had a problem\n');
    j
    size(File.Distance)
    Indices(j)
    a
end

  end
  NewIndices = [NewIndices find(a)];
end

NewIndices = NewIndices(NewIndices > 0);
NewIndices = NewIndices(NewIndices <= File.NumNT);

NewIndices = unique(NewIndices);
NewIndices = sort(NewIndices);
NewIndices = NewIndices(NewIndices>0); % Anton 1/8/2011 Prevents overflow problem
NewIndices = NewIndices(NewIndices<=File.NumNT); % Anton 1/14/2011 Prevents overflow problem