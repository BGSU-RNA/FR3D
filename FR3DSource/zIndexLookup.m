% zIndexLookup(File,Num,Chain) finds the base index for base having nucleotide 
% number Num and, if specified, Chain
% Num can be a cell array of N nucleotide numbers
% Any entry of that cell array can use the notation '1830:1835' for a range
% Any entry can also use commas to separate nucleotide numbers
% Nucleotide numbers can be followed by (A) or _A to indicate chain A
% Someday: Base letters can be specified to help identify the correct chain

% ind is a Cx1 vector of the first match to each, the easy answer
% allchains is a Cx1 cell array.  Each cell lists all chains that match the
% corresponding nucleotide number

function [ind,allchains] = zIndexLookup(File,Num,Chain,Verbose)

if nargin < 3,
  for k = 1:length(Num),
    Chain{k} = '';
  end
end

if nargin < 4,
  Verbose = 1;
end

if strcmp(class(Num),'char'),
  Num = {Num};
end

if strcmp(class(Chain),'char'),
  Chain = {Chain};
end

% split multiple specifications into separate cells

t = 1;
for k = 1:length(Num),
  N = regexprep(Num{k},'\t',' ');
  N = regexprep(N,'([ACGU])\s*([123456789])','$1$2'); % join base & number
  N = regexprep(N,';| ',',');         % replace other delimeters with commas
  while strfind(N,',,'),              % replace double commas with single
    N = regexprep(N,',,',',');
  end
  a = [1 1+strfind(N,',')];               % locations of commas
  for i = 1:(length(a)-1),
    Numb{t} = N(a(i):(a(i+1)-2));
    Chai{t} = Chain{k};
    t = t + 1;
  end
  Numb{t} = N(a(end):end);
  Chai{t} = Chain{k};
  t = t + 1;
end

% check for chain indicated in parentheses

for k = 1:length(Numb),
  if Numb{k}(1) == '(' || Numb{k}(1) == '_',         % chain only
    Chai{k} = Numb{k}(2);                            % extract the chain
    Numb{k} = '#';                                   % special code
  elseif ~isempty(strfind(Numb{k},'(')),
    a = strfind(Numb{k},'(');
    b = strfind(Numb{k},')');
    Chai{k} = Numb{k}(a(1)+1:b(1)-1);     % extract chain
    if length(a) == 1,                    % one chain specified      
      Numb{k} = [Numb{k}(1:a(1)-1) Numb{k}(b(1)+1:end)];
    elseif length(a) == 2,                % range with two chains specified
      Numb{k} = [Numb{k}(1:a(1)-1) Numb{k}(b(1)+1:a(2)-1)];
    end
  end
end

% check for element or chain indicated by underscore

for k = 1:length(Numb),
  if any(Numb{k}(1) == 'bdefhijklmnopqrstvwxyz'),
    % this is an element, leave it alone
  else                      % check for chain indicated by underscore
    if ~isempty(strfind(Numb{k},'_')),
      a = strfind(Numb{k},'_');
      Chai{k} = Numb{k}(a(end)+1:end);        % extract chain
      Numb{k} = Numb{k}(1:a(end)-1);          % remove chain reference
    end
  end
end

ind = [];
  
% if File is a text string (filename), load the file

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

Numbers = cat(1,{File.NT(:).Number});

allchains = [];

for k = 1:length(Numb)                      % loop through nucleotide numbers
  if ~isempty(strfind(Numb{k},':')),
                                            % a range of numbers
    n = Numb{k};                            % kth specified number or range
    i = strfind(n,':');                     % find the colon
    Numb1 = n(1:(i-1));                     % first nucleotide number
    p = LookUpOne(File,Numbers,Numb1,Chai{k},Verbose);
    Numb2 = n((i+1):length(n));             % second nucleotide number
    q = LookUpOne(File,Numbers,Numb2,Chai{k},Verbose);

    m = length(ind);

    if length(p) > 0 && length(q) > 0,      % both numbers are valid
      if p(1) < q(1),
        ind = [ind p(1):q(1)];
      else
        ind = [ind p(1):-1:q(1)];
      end
    end

    mm = length(ind);
    ch = [];
    c  = 1;

    for j = 1:length(p),
      for jj = 1:length(q),
        if File.NT(p(j)).Chain == File.NT(q(jj)).Chain,
          ch{c} = File.NT(p(j)).Chain;
          c = c + 1;
        end
      end
    end

    for j = (m+1):mm,
      allchains{j} = ch;
    end
  elseif Numb{k}(1) == '#',
    ind = find(cat(2,File.NT.Chain) == Chai{k});
  elseif strcmpi(Numb{k},'all'),
    ind = 1:length(File.NT);                  % all nucleotides
  elseif any(Numb{k}(1) == 'bdefhijklmnopqrstvwxyz') && isfield(File.NT(1),'Hierarchy'),
    c = 1;                        % counter for number of indices found
    newind = [];
    for n = 1:length(File.NT),
      p = find(ismember(File.NT(n).Hierarchy,Numb{k}));
      if ~isempty(p),
        newind(c) = n;
        c = c + 1;
      end
    end

    if ~isempty(newind),
      for j = 1:length(newind),
        allchains{length(ind)+j} = File.NT(newind(j)).Chain;
      end
      ind = [ind newind];
    end
  else
    L = LookUpOne(File,Numbers,Numb{k},Chai{k},Verbose);

    if length(L) > 0,
      ind = [ind L(1)];                       % add the first hit to the list

      m = length(ind);
      ch = [];

      for i=1:length(L),
        ch{i} = File.NT(L(i)).Chain;
      end

      allchains{m} = ch;
    end
  end
end

return

% check for indicated base - add to LookUpOne in the future

for k = 1:length(Numb),
  Numb{k} = upper(Numb{k});
  aa = strfind(Numb{k},'A');
  ac = strfind(Numb{k},'C');
  ag = strfind(Numb{k},'G');
  au = strfind(Numb{k},'U');
  Base{k} = '';
  if ~isempty(aa),
    Numb{k} = strrep(Numb{k},'A','');
    Base{k} = 'A';
  end
  if ~isempty(ac),
    Numb{k} = strrep(Numb{k},'C','');
    Base{k} = 'C';
  end
  if ~isempty(ag),
    Numb{k} = strrep(Numb{k},'G','');
    Base{k} = 'G';
  end
  if ~isempty(au),
    Numb{k} = strrep(Numb{k},'U','');
    Base{k} = 'U';
  end
end

%-------------------------------------------------------------------------
function [ind] = LookUpOne(File,Numbers,N,Chain,Verbose)

    if any(N(1) == 'ACGU'),
      N = N(2:end);
    end

    % later, add the ability to infer the chain from the base specified

    ind = [];
    p = find(ismember(Numbers,N));
    if length(p) == 0,
      if Verbose > 0,
        fprintf('Could not find nucleotide %s in %s\n',N,File.Filename);
      end
    elseif length(p) == 1 & length(Chain) == 0, % one match, no chain specified
      ind = [ind p];
    elseif length(p) > 1 & length(Chain) == 0,% two matches, no chain specified
      ind = [ind p];
      if Verbose > 0,
        fprintf('Multiple matches found for %s in %s, consider specifying a chain\n', N, File.Filename);
      end
      for a = 1:length(ind),
        if Verbose > 0,
          fprintf('Nucleotide %s%s Chain %5s Index %5d\n', File.NT(ind(a)).Base, File.NT(ind(a)).Number, File.NT(ind(a)).Chain, ind(a));
        end
      end
    elseif length(Chain) > 0,                    % chain specified
      c = 0;
      for j = 1:length(p),
        if strcmp(File.NT(p(j)).Chain,Chain),
          ind = [ind p(j)];
          c = c + 1;
        end
      end
      if c == 0,
        if Verbose > 0,
          fprintf('Could not find nucleotide %s in chain %s in %s\n',N,Chain,File.Filename);
        end
      elseif c > 1,
        if Verbose > 0,
          fprintf('Multiple matches found for %s in chain %s in %s\n', N,Chain,File.Filename);
        end
      end
    end
