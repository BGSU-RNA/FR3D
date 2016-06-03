% xBorderSS identifies nucleotides which border a single-stranded region by making nested canonical cWW basepairs

function [File] = xBorderSS(File,Verbose)

if nargin < 2,
  Verbose = 0;
end

if File.NumNT > 0,
  E = triu(fix(abs(File.Edge))==1);     % cWW pairs
  [i,j] = find(E);                      % indices of NT's making cWWs

  % Paircode list
  % 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC
  % 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

  %  paircode = 4*(N2.Code-1) + N1.Code;           % AA is 1, CA is 2, etc.

  Nested = zeros(1,length(i));          % place to stored nested canonical cWWs
  for u = 1:length(i),
    pc = 4*(File.NT(j(u)).Code-1) + File.NT(i(u)).Code;
    if File.Crossing(i(u),j(u)) == 0 && any(pc == [4 7 10 12 13 15]), % canonical
      Nested(u) = 1;
    end
  end

  k = find(Nested);
  i = i(k);                             % indices making nested cWWs
  j = j(k);

  cWWPartner = zeros(1,File.NumNT);     % non-zero if there is a nested canonical cWW
  cWWPartner(i) = j;
  cWWPartner(j) = i;

  a = [i; j];                           % all indices of nested canonical cWWs
  a = sort(a);                          % put them into one list

  H = sparse([],[],[],File.NumNT,File.NumNT);    % matrix to indicate borderSS
  bSSPartnerBefore = zeros(1,File.NumNT);     % non-zero for borderSS partners
  bSSPartnerAfter  = zeros(1,File.NumNT);     % non-zero for borderSS partners

  H(1,min(a)) = 1;                      % 5' end of the chain

  for k = 1:(length(a)-1),              % run through all NT making nested cWW
    if a(k+1) - a(k) > 1,               % where there is a gap in the number,
      if strcmp(File.NT(a(k)).Chain, File.NT(a(k+1)).Chain),  % and same chain,
        H(a(k),a(k+1)) = 1;               % these two border a single-stranded region
        H(a(k+1),a(k)) = 1;               % these two border a single-stranded region
      end
    end
  end

  % Record border SS partners; each nucleotide can have one before and/or after in the chain

  [i,j] = find(triu(H));
  bSSPartnerAfter(i) = j;
  bSSPartnerBefore(j) = i;

  % Identify self-complementary internal loops and remove them from borderSS

  p = find(bSSPartnerAfter);              % indices with non-zero values

  for k = 1:length(p),
    a = p(k);                             % first nucleotide of a pair
    d = cWWPartner(a);                    % cWW partner of a
    if d > a,                             % no need to check a loop multiple times
      b = bSSPartnerAfter(a);             % borderSS partner of a further down the chain
      if b ~= d,                          % this is not an HL
        c = cWWPartner(b);                % nucleotide on other strand
        if c == bSSPartnerBefore(d),      % this is an IL
          if b-a == d-c,                  % symmetric IL
%[File.NT(a).ID ' ' File.NT(b).ID ' ' File.NT(c).ID ' ' File.NT(d).ID]
            j = [a:b c:d];                % all indices of nts in this loop
            e = fix(abs(File.Edge));      % all basepair interactions
            e(a,d) = 0;                   % ignore flanking cWW interactions
            e(b,c) = 0;                   % ignore flanking cWW interactions
            e = fix(abs(File.Edge(j,j))); % all basepair interactions in this loop

            if sum(sum((e > 1) .* (e < 13))) == 0, % no non-cWW basepairs here
              selfcomp = 1;               % assume self complementary
              for m = 1:(b-a-1),          % loop through nucleotides
                pc = 4*(File.NT(a+m).Code-1) + File.NT(d-m).Code;
                if ~any(pc == [4 7 10 12 13 15]),
                  selfcomp = 0;           % non-complementary pair found
                end
              end
              if selfcomp == 1,           % all self complementary
                if Verbose > 0,
                  zShowInteractionTable(File,j);
                end
                H(a,b) = 0;               % remove borderSS relations for this loop
                H(b,a) = 0;
                H(c,d) = 0;
                H(d,c) = 0;
              end
            end
          end
        end
      end
    end
  end

  % Record border SS partners; each nucleotide can have one before and/or after in the chain

  [i,j] = find(triu(H));
  bSSPartnerAfter(i) = j;
  bSSPartnerBefore(j) = i;

  % Find adjacent bases making cWWs with BorderSS NTs

  p = find(bSSPartnerBefore);             % indices at the end of a single-stranded region

  for k = 1:length(p),
    b = p(k);                             % nucleotide at the end of a single-stranded region
    bb = b;                               % starting nucleotide
    c = cWWPartner(b);                    % basepairing partner of b
    if c+1 <= File.NumNT,                 % make sure we do not go past the end of the file
      d = cWWPartner(c+1);                  % cWW partner of the next nucleotide after c, if any
      while d > 0 && d ~= bb-1,              % if c+1 makes a nested cWW and we have not closed the loop,
        H(c,c+1) = 1;
        H(c+1,c) = 1;
        if Verbose > 0,
          fprintf('Adding %-16s - %-16s to BorderSS list, cWWs with %-16s and %-16s (difference of %5d)\n',File.NT(c).ID, File.NT(c+1).ID, File.NT(b).ID, File.NT(d).ID,b-d);
        end
        b = c+1;
        c = d;
        if c+1 <= File.NumNT,
          d = cWWPartner(c+1);
        else
          d = 0;
        end
      end
    end
  end

  % Record border SS partners; each nucleotide can have one before and/or after in the chain

  [i,j] = find(triu(H));
  bSSPartnerAfter(i) = j;
  bSSPartnerBefore(j) = i;

  % Add 5' ends of each chain to BorderSS list
  CurrentChain = '';

  for n = 1:(File.NumNT-1),
    if ~strcmp(File.NT(n).Chain,CurrentChain),
      a = n;                                    % a is the first nucleotide of a chain

      CurrentChain = File.NT(n).Chain;

      while cWWPartner(a) == 0 && a < File.NumNT && strcmp(File.NT(a+1).Chain,File.NT(n).Chain),
        a = a + 1;
      end

      if a > n,
        H(a,n) = 1;
        H(n,a) = 1;

        if Verbose > 0,
          fprintf('Adding %-16s - %-16s to BorderSS list because it is at the 5'' end of chain %s\n',File.NT(n).ID, File.NT(a).ID,CurrentChain);
        end
      else
        H(a,a) = 1;
        if Verbose > 0,
          fprintf('Adding %-16s - %-16s to BorderSS list because it is at the 5'' end of chain %s and makes a nested cWW\n',File.NT(a).ID, File.NT(a).ID,CurrentChain);
        end
      end
    end
  end

  % Add 3' ends of each chain to BorderSS list
  CurrentChain = '';

  for n = 2:File.NumNT,
    if ~strcmp(File.NT(n).Chain,CurrentChain) || n == File.NumNT,
      if n == File.NumNT,
        a = n;                                    % a is the last nucleotide in the file
      else
        a = n-1;                                  % a is the last nucleotide of a chain
      end
      aa = a;                                     % original value of a
      while cWWPartner(a) == 0 && a > 1 && strcmp(File.NT(a-1).Chain,File.NT(aa).Chain),
        a = a - 1;
      end
      if a < aa,
        H(a,aa) = 1;
        H(aa,a) = 1;
        if Verbose > 0,
          fprintf('Adding %-16s - %-16s to BorderSS list because it is at the 3'' end of chain %s\n',File.NT(a).ID, File.NT(aa).ID,CurrentChain);
        end
      elseif a > 1,
        H(a,a) = 1;
        if Verbose > 0,
          fprintf('Adding %-16s - %-16s to BorderSS list because it is at the 3'' end of chain %s and makes a nested cWW\n',File.NT(a).ID, File.NT(a).ID, File.NT(a).Chain);
        end
      end

      CurrentChain = File.NT(n).Chain;

    end
  end

  % Remove BorderSS that crosses two chains or two symmetry operators

  [i,j] = find(triu(H));
  for k = 1:length(i),
    a = i(k);
    b = j(k);
    if ~strcmp(File.NT(a).Chain,File.NT(b).Chain),
      H(a,b) = 0;
      H(b,a) = 0;
      if Verbose > 0,
        fprintf('Removed BorderSS between %s and %s because they are in different chains\n',File.NT(a).ID,File.NT(b).ID);
      end
    end

    afields = zNTIDFields(File.NT(a).ID);
    bfields = zNTIDFields(File.NT(b).ID);

    if ~strcmp(afields{9},bfields{9}),
      H(a,b) = 0;
      H(b,a) = 0;
      if Verbose > 0,
        fprintf('Removed BorderSS between %s and %s because they have different symmetry operators\n',File.NT(a).ID,File.NT(b).ID);
      end
    end

  end

  File.Flank = H;                       % store


  % Iterate through hairpins, internal loops, and junction loops

  if Verbose > 1,
    % Record border SS partners; each nucleotide can have one before and/or after in the chain

    [i,j] = find(triu(H));
    bSSPartnerAfter(i) = j;
    bSSPartnerBefore(j) = i;

    Visited(i) = 0;
    Visited(j) = 0;

    i = sort(i);

    for k = 1:length(i),
      aa = i(k);
      a = i(k);
      c = 0;
      if Visited(a) == 0,
        while a > 0,
          fprintf(' %s', File.NT(a).ID);
          Visited(a) = 1;
          b = bSSPartnerAfter(a);
          if b > 0,
            Visited(b) = 1;
            fprintf(' bSS %s',File.NT(b).ID);
            c = cWWPartner(b);
            if c > 0,
              Visited(c) = 1;
              fprintf(' cWW %s',File.NT(c).ID);
              a = bSSPartnerAfter(c);
            end
          end
        end
      end

      fprintf('\n');

    end
  end

else

  File.Flank = [];

end

