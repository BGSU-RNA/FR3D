
function [void] = zUpdatePDBInfo(File,Verbose)

load PDBInfo               % load current data file of information on each PDB file

i = find(ismember(t(:,1),File.Filename));
if length(i) > 1,
  fprintf('zUpdatePDBInfo:  more than one entry for %s\n',File.Filename);
  i = i(1);
end

if isempty(i),
  i = size(t,1) + 1;                % add this file at the end
end

if ~isempty(File.NT),                      % if it has nucleotides, update the data and save
  t{i,1} = File.Filename;
  if isfield(File.Info,'Descriptor'),
    t{i,2} = File.Info.Descriptor;
  end
  if isfield(File.Info,'ExpTechnique'),
    t{i,3} = File.Info.ExpTechnique;
  end
  if isfield(File.Info,'ReleaseDate'),
    t{i,4} = File.Info.ReleaseDate;
  end
  if isfield(File.Info,'Author'),
    t{i,5} = File.Info.Author;
  end
  if isfield(File.Info,'Keywords'),
    t{i,6} = File.Info.Keywords;
  end
  if isfield(File.Info,'Resolution'),
    t{i,7} = File.Info.Resolution;
  end
  if isfield(File.Info,'Source'),
    t{i,8} = File.Info.Source;
  end
%  t{i,9} is not used
%  t{i,10} is set by zFileRedundancy to indicate the representative of the class at all resolutions

  LC = File.LongestChain;

  t{i,11} = cat(2,File.NT(LC(1):LC(2)).Base);    % bases in longest chain
  t{i,12} = File.BestChains;                     % characters of the best chain(s)
  t{i,13} = File.NT(LC(1)).Chain;                % identifier of the longest chain

File.NT(LC(1)).Chain



  n(i,1) = File.Info.Resolution;
  n(i,2) = length(File.NT);                % store the number of NT

  E  = abs(triu(File.Edge));
  n(i,3) = full(sum(sum((E > 0) .* (E < 13)))); % number of pairs

  n(i,4) = length(t{i,11});         % number of nucleotides in longest chain
  n(i,5) = LC(1);                   % starting index of longest chain
  n(i,6) = LC(2);                   % end index of longest chain

% n(i,7)    Number of nucleotides in NR chains
% n(i,8)    Number of basepairs in NR chains
% n(i,9)    Number of cWW basepairs in structure
% n(i,10)   Number of non-cWW basepairs in structure

  save(which(['FR3DSource' filesep 'PDBInfo.mat']),'n','t'); % Matlab version 7
  fprintf('zUpdatePDBInfo:  Updated PDBInfo.mat\n');
end




return

current = 1;
load PDBInfo
for i = current:length(t(:,1)),

  current = i;

  File = zAddNTData(t{i,1},0,[],1);          % load file
  fprintf('Updating structure %d, which is %s and has %d nucleotides\n', i, t{i,1}, length(File.NT));

  zUpdatePDBInfo(File);

end

return

load PDBInfo
current = size(t,1);
for i = current:-1:1,

  current = i;

  File = zAddNTData(t{i,1},0,[],1);          % load file
  fprintf('Updating structure %d, which is %s and has %d nucleotides\n', i, t{i,1}, length(File.NT));

  zUpdatePDBInfo(File);

end



return

for i = size(t,1):-1:1,
  if isempty(t{i,13}),
    i
    File = zAddNTData(t{i,1},0,[],1);          % load file
    File
    pause

    if File.NumNT <= 0 && ~isempty(File.LongestChain),
      File.LongestChain = [];
      zSaveNTData(File);
    elseif File.NumNT > 0,
      File
      pause
    end
  end
end

return


for i = current:length(t(:,1)),

  fprintf('Updating structure %d, which is %s\n', i, t{i,1});

  current = i;

  try
    File = zAddNTData(t{i,1},0,[],1);          % load file
  catch
    delete(['PrecomputedData' filesep t{i,1} '.mat']);
    File = zAddNTData(t{i,1},0,[],1);          % load file
  end

  if length(File.NT) > 1,
    c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
    File.Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
    d = sort(nonzeros(File.Distance));

    if length(d) > 1,
      if d(min(10,length(d))) < 1,
        File = zAddNTData(t{i,1},4,[],1);     % fix NMR files
      end
    end
  end

  if ~isempty(File.NT),                      % if it has nucleotides,
    n(i,2) = length(File.NT);                % store the number of NT

    E  = abs(triu(File.Edge));
    n(i,3) = full(sum(sum((E > 0) .* (E < 13)))); % number of pairs

    LC = File.LongestChain;

    t{i,11} = cat(2,File.NT(LC(1):LC(2)).Base);    % bases in longest chain
    t{i,12} = File.BestChains;        % characters of the best chain(s)

    n(i,4) = length(t{i,11});         % number of nucleotides in longest chain
    n(i,5) = LC(1);                   % starting index of longest chain
    n(i,6) = LC(2);                   % end index of longest chain

    if Verbose > 1,
      fprintf('All      %s\n',cat(2,File.NT.Base));
      fprintf('Longest  %s\n',t{i,11});
    end
  end

%  zSaveNTData(File);

  if mod(i,100) == 0,
    save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
  end

end
