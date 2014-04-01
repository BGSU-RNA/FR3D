% xNeighborhood(File,Indices,v,MaxDiff,MaxInsert) returns Indices plus indices "near"
% the given indices, in a way determined by the vector v

function [NewIndices] = xNeighborhood(File,Indices,v,strandnumber)

N = length(Indices);                                 % number of nucleotides

NewIndices = [];                                     % always start here

if max(v) > 0,
  if v(1) > 0,                                            % fill in strands
    if nargin < 4,
      for a = 1:N,
        for b = (a+1):N,
          if abs(double(Indices(a))-double(Indices(b))) <= 5,
            NewIndices = [NewIndices min(Indices(a),Indices(b)):max(Indices(a),Indices(b))];
          end
        end
      end
    else
      for a = 1:N,
        for b = (a+1):N,
          if strandnumber(a) == strandnumber(b),
            NewIndices = [NewIndices min(Indices(a),Indices(b)):max(Indices(a),Indices(b))];
          end
        end
      end
    end
    NewIndices = unique(NewIndices);
    NewIndices = setdiff(NewIndices,Indices);        % keep distinct from original indices
  end

  if v(2) > 0,                                            % extend strands 1
    NewIndices = [Indices-1 Indices+1];
    NewIndices = NewIndices(NewIndices > 0);
    NewIndices = NewIndices(NewIndices <= length(File.NT));
  end

  if v(3) > 0,                                            % extend strands 2
    NewIndices = [Indices-1 Indices+1 Indices-2 Indices+2];
    NewIndices = NewIndices(NewIndices > 0);
    NewIndices = NewIndices(NewIndices <= length(File.NT));
  end

  if v(4) > 0,                                            % extend strands 3
    NewIndices = [Indices-1 Indices+1 Indices-2 Indices+2 Indices-3 Indices+3];
    NewIndices = NewIndices(NewIndices > 0);
    NewIndices = NewIndices(NewIndices <= length(File.NT));
  end

  if v(5) > 0,                                          % add basepairs
    for n = 1:N,
      e = abs(File.Edge(Indices(n),:));              % interactions with n
      j = find( (e > 0) .* (e < 14) );               % basepairing with n
      NewIndices = [NewIndices j];
    end
    for n = 1:length(NewIndices),                    % add what NewIndices pair with
      e = abs(File.Edge(NewIndices(n),:));           % interactions with n
      j = find( (e > 0) .* (e < 14) );               % basepairing with n
      NewIndices = [NewIndices j];
    end
    NewIndices = unique(NewIndices);
    NewIndices = setdiff(NewIndices,Indices);        % keep distinct from original indices
  end

  if v(6) > 0,                                          % add stacks
    for n = 1:length(NewIndices),                    % add what NewIndices stack on
      e = abs(File.Edge(NewIndices(n),:));           % interactions with n
      j = find( (e > 20) .* (e < 24) );               % stacking with n
      NewIndices = [NewIndices j];
    end
    for n = 1:N,
      e = abs(File.Edge(Indices(n),:));              % interactions with n
      j = find( (e > 20) .* (e < 24) );               % stacking with n
      NewIndices = [NewIndices j];
    end
    NewIndices = unique(NewIndices);
    NewIndices = setdiff(NewIndices,Indices);        % keep distinct from original indices
  end

  if v(7) > 0,                                          % add stacks
    for n = 1:length(Indices),                          % add base-backbone interactions with Indices
      j1 = find(File.BasePhosphate(Indices(n),:) > 0);
      j2 = find(File.BasePhosphate(:,Indices(n)) > 0);
      j3 = find(File.BaseRibose(Indices(n),:) > 0);
      j4 = find(File.BaseRibose(:,Indices(n)) > 0);
      NewIndices = [NewIndices j1 j2' j3 j4'];
    end
    NewIndices = unique(NewIndices);
    NewIndices = setdiff(NewIndices,Indices);        % keep distinct from original indices
  end

  if v(8) > 0,
    if ~isfield(File,'Distance'),
      File.Distance = [];
    end

    [s,t] = size(File.Distance);

    if isempty(File.Distance) || s < File.NumNT,
      c = 9999999 * ones(File.NumNT,3);
      for i = 1:length(File.NT),
        if ~isempty(File.NT(i).Center),
          c(i,:) = File.NT(i).Center;
        end
      end

      File.Distance = sparse([],[],[],File.NumNT,File.NumNT);
      File.Distance = zMutualDistance(c,v(8));           % compute distances needed
    end

    a = zeros(1,File.NumNT);
    for j=1:length(Indices),
      try
          a = a + (File.Distance(Indices(j),:) < v(8)) .* (File.Distance(Indices(j),:) > 0);
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

  NewIndices = NewIndices(NewIndices > 0);           % Anton 1/8/2011 Prevents overflow problem
  NewIndices = NewIndices(NewIndices <= File.NumNT); % Anton 1/14/2011 Prevents overflow problem

  NewIndices = unique(NewIndices);
  NewIndices = sort(NewIndices);
  NewIndices = setdiff(NewIndices,Indices);        % keep distinct from original indices

  removeindices = [];
  for j = 1:length(NewIndices),
    if isempty(File.NT(NewIndices(j)).Fit),
      removeindices = [removeindices NewIndices(j)];
    end
  end
  NewIndices = setdiff(NewIndices,removeindices);

end