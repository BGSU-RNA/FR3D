
function [File] = zBackboneConformation(File,Verbose)

if nargin < 2,
  Verbose = 0;
end

File.Backbone = sparse([],[],[],length(File.NT),length(File.NT));

if exist('chiropraxis.jar') == 2 ...       % if chiropraxis is here
   && File.NumNT > 1 ...                   % and there is more than one NT
   && exist(File.PDBFilename) == 2,        % and the PDB file is available

  d = which(File.PDBFilename);
  e = which('chiropraxis.jar');
  g = which('suitename.0.3.070628.win.exe');

  if Verbose > 0,
    fprintf('Using dangle and suitename to classify backbone conformations.\n');
  end

  system(['java -cp "' e '" chiropraxis.dangle.Dangle rnabb < "' d '" > ' File.Filename '_dangle.txt']);

  system(['"' g '" < ' File.Filename '_dangle.txt > ' File.Filename '_suitename.txt']);

%  fid = fopen('tempbat.bat','w');
%  fprintf(fid,'echo off\n');
%  fprintf(fid,'java -cp "%s" chiropraxis.dangle.Dangle rnabb < "%s" > %s_dangle.txt\n', e, d, File.Filename);
%  fprintf(fid,'"%s" < %s_dangle.txt > %s_suitename.txt\n', g, File.Filename, File.Filename);
%  fprintf(fid,'pause\n');
%  fclose(fid);
%  !tempbat.bat

  % ------------------------------------------- read data written by suitename

  fid = fopen([File.Filename '_suitename.txt'],'r');

  if fid > 0

  clear T
  T = [];
  L = 1;
  c = 1;

  while L > -1
    L = fgets(fid);
    if L > -1
      if L(1) == ':',
        T{c} = L;
        c = c + 1;
      else
        L = -1;                                   % stop reading file
      end
    end
  end

  else

  fprintf('zBackboneConformation could not open file %s\n', [File.Filename '_suitename.txt']);

  end

  fclose(fid);

  %123456789012345678901234567890
  %:1:A: 388: :  G 33 p 1a 0.630

  File.Backbone = sparse([],[],[],length(File.NT),length(File.NT));

  zBackboneCodes;

  Numbers = cat(1,{File.NT(:).Number});

  c = 0;                                      % current index in the file

  

  for t = 1:length(T),
  a = T{t};

  h = strfind(a,':');

  if a(h(4)+1) ~= ' ',
    b = [a((h(3)+1):(h(3)+4)) a(h(4)+1)];
  else
    b = a((h(3)+1):(h(3)+4));
  end

  i = LookUpOne(File,Numbers,strrep(b,' ',''),a(h(2)+1),0);

  if length(i) > 1,
    w = find(i > c);
    i = i(w(1));
  end

  if ~isempty(i),
    c = i;
  end

  if isempty(i),
    if Verbose > 0,
      fprintf('%s has no corresponding nucleotide in FR3D\n', a(1:29));
    end
  elseif Verbose > 2,
    fprintf('%s becomes %s%s\n', a(1:29), File.NT(i).Base, File.NT(i).Number);
  end

  bbc = find(ismember(Codes,a((h(5)+10):(h(5)+11))));   % look up FR3D's internal code

  if isempty(bbc),
    fprintf('Found a new code, %s\n', a((h(5)+10):(h(5)+11)));
  elseif bbc < length(Codes) && ~isempty(i),
    if i > 1,
      if File.Covalent(i-1,i) == 1,             % linked to previous NT
        File.Backbone(i-1,i) = bbc;
      end
    end
  end
  end


%  delete('tempbat.bat');

  delete([File.Filename '_dangle.txt']);
  delete([File.Filename '_suitename.txt']);

else
%  fprintf('Cannot determine backbone conformations for %s\n', File.Filename);
end

return

% determine necessary codes

for t = 1:length(T),
  g{t} = T{t}(22:23);
end
G = unique(g);


% now read the successful output of suitename, convert the codes to numbers, store them, match the codes to numbers in zEdgeText, and implement this in xPairwiseScreen.  Also make it possible to see the backbone conformations for candidates.  Work, but not so bad.


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
