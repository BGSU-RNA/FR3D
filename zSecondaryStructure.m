% zSecondaryStructure(File,NTNumber) displays the secondary structure of a
% molecule as reflected by the interaction matrix, starting at NTNumber and
% ending with whather base NTNumber interacts with

function [File,Interact] = zSecondaryStructure(File,NTNumber,LastNTNumber,Interact)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = {File};
  File = zGetNTData(Filename,0);
end

N = length(File.NT);                           % number of nucleotides

% if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char'),
  NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell'),
  FirstIndex = zIndexLookup(File,NTNumber);
else
  FirstIndex = NTNumber;
end

if nargin > 2,
  if strcmp(class(LastNTNumber),'char'),
    LastNTNumber = {LastNTNumber};
  end
  if strcmp(class(LastNTNumber),'cell'),
    LastIndex = zIndexLookup(File,LastNTNumber);
  else
    LastIndex = LastNTNumber;
  end
else
  LastIndex = N;
end

E = abs(fix(File.Edge));                   % don't distinguish subcategories

if nargin < 4,
 Basepairs = (E < 15) .* (E ~= 0);         % basepairing only
 for a = 1:N,                              % loop through nucleotides
  k = find(Basepairs(a,:));                % find indices of interacting bases
  [y,L] = sort(E(a,k));                    % sort by interaction category
  Interact{a}.Categ = E(a,k(L));           % store categories
  Interact{a}.Index = k(L);                % store indices of interacting bases
 end
end

a = FirstIndex(1);                         % first index

if length(Interact{a}.Categ) > 0,          % if a interacts with something,
  b = Interact{a}.Index(1);                % get index
  c = Interact{a}.Categ(1);                % and category of interaction
  if (b < a) & (c == 1),                   % canonical pair with lower index
%    a = b;                                 % swap a and b
  end
end

b = Inf;                                   % index of base it interacts with
c = 2;                                     % category of interaction

A  = a;
B  = LastIndex;
AA = a;
BB = a;

while (a < B) & (a <= N), % while not the end of the loop,

  if length(Interact{a}.Categ) > 0,        % if a interacts with something,
    b = Interact{a}.Index(1);              % get index 
    c = Interact{a}.Categ(1);              % and category of interaction
    t = zEdgeText(File.Edge(a,b));         % edge notation
    if (a < b),                            % if b comes after a

      % ---------------------------------- Check for junction

      if (sum(sum(E((b+1):(BB-1),(b+1):(BB-1)) == 1)) > 0) & ...
         (sum(sum(E((a+1):(b-1),(a+1):(b-1))   == 1)) > 0),
        fprintf('\nJunction\n\n');
        fprintf('Loop 1 - Nucleotides %s to %s, length %3d\n',File.NT(a).Number,File.NT(b).Number,b+1-a);
        zSecondaryStructure(File,a,b,Interact);
        fprintf('Loop 2 - Nucleotides %s to %s, length %3d\n',File.NT(b+1).Number,File.NT(BB-1).Number,BB-b);
        zSecondaryStructure(File,b+1,BB-1,Interact);

        return
      end

      % ---------------------------------- Check for insertions on right

      if (b < B-1),
        if sum(sum(E((b+1):(B-1),(b+1):(B-1)) == 1)) == 0,
          for e = (B-1):-1:(b+1),
            fprintf('           %1s%4s', File.NT(e).Base, File.NT(e).Number);
            if length(Interact{e}.Categ) > 0,
             fprintf('   Right ');
             for k=1:length(Interact{e}.Categ),
              bb = Interact{e}.Index(k);
              cc = Interact{e}.Categ(k);
              te = zEdgeText(File.Edge(e,bb));
              fprintf('%1s%4s',File.NT(e).Base,File.NT(e).Number);
              fprintf(' %4s %1s%4s | ',te,File.NT(bb).Base,File.NT(bb).Number);
             end
            end
            fprintf('\n');
          end
        end
      end

      % ---------------------------------- Check if b is out of sequence

      if ((sum(sum(E((a+1):(b-1),(a+1):(b-1)) == 1)) == 0) & ...
          (sum(sum(E((a+1):B,    (a+1):B)     == 1))  > 0)) | ...
         (b > B),
        fprintf('%1s%4s',File.NT(a).Base,File.NT(a).Number); % display a
        fprintf('           ');       % move to next column
        d = 1;
      else
        % ---------------------------------- Show interaction between a and b

        fprintf('%1s%4s',File.NT(a).Base,File.NT(a).Number); % display a
        fprintf(' %4s %1s%4s', t, File.NT(b).Base, File.NT(b).Number);
        d = 2;
        A = a;                               % a of last pair
        B = b;                               % b of last pair
        if c == 1,
          AA = a;                            % a of last class 1 pair
          BB = b;                            % b of last class 1 pair
        end
      end
    else
      fprintf('%1s%4s',File.NT(a).Base,File.NT(a).Number); % display a
      fprintf('           ');       % move to next column
      d = 1;
    end
  else
    fprintf('%1s%4s',File.NT(a).Base,File.NT(a).Number); % display a
    d = 0;
  end

  % ---------------------------------- Show additional interactions a has

  if (((d == 2) & (length(Interact{a}.Categ) > 1)) | (d == 1)),
    fprintf('   Left  ');
    for k = d:length(Interact{a}.Categ),
      bb = Interact{a}.Index(k);
      cc = Interact{a}.Categ(k);
      te = zEdgeText(File.Edge(a,bb));
      fprintf('%1s%4s',File.NT(a).Base,File.NT(a).Number);
      fprintf(' %4s %1s%4s | ',te,File.NT(bb).Base, File.NT(bb).Number);
    end
  end

  % ---------------------------------- Show additional interactions b has

  if d == 2,                              % if there was a primary interaction
    if length(Interact{b}.Categ) > 1,     % and b makes more than one interact
      fprintf('   Right ');
      for k=1:length(Interact{b}.Categ),
        bb = Interact{b}.Index(k);
        cc = Interact{b}.Categ(k);
        te = zEdgeText(File.Edge(b,bb));
        if bb ~= a,
          fprintf('%1s%4s',File.NT(b).Base,File.NT(b).Number);
          fprintf(' %4s %1s%4s | ',te,File.NT(bb).Base,File.NT(bb).Number);
        end
      end
    end
  end

  fprintf('\n');
  a = a + 1;
end

fprintf('\n');
