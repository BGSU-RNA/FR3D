% zEdgeMakesMultiplePairs identifies nucleotides which have been classified as using the same edge in more than one basepair and chooses the best basepair among them

function [File] = zEdgeMakesMultiplePairs(File,Verbose)

% File = zAddNTData('2avy');
% File = zAddNTData('Nonredundant_2009-05-14_list');

if nargin < 2,
  Verbose = 1;                  % list the offending pairs
  Verbose = 3;                  % show the nucleotides and wait for key press
  Verbose = 0;                  % don't say anything
  Verbose = 2;                  % notify of the change
end

load PairExemplars

rr = 0;                             % total ncWW removed
nn = 0;                             % total original ncWW

for f = 1:length(File),

  % ----------------------------- remove ncWW if also making cWW

  E = abs(fix(File(f).Edge));
  [i,j] = find(triu(E == 101));     % ncWW

  r = 0;                            % number removed from File(f)
  nn = nn + length(i);              % total number of cWW found

  for k = 1:length(i),              % loop through ncWW
    if any(E(i(k),:) == 1) && any(E(j(k),:) == 1),  % both make cWW
      File(f).Edge(i(k),j(k)) = 0;  % remove ncWW
      File(f).Edge(j(k),i(k)) = 0;  % remove ncWW

      r = r + 1;
      rr = rr + 1;

      if Verbose == 3,
        ii = [max(i(k)-1,1) i(k) min(i(k)+1,File(f).NumNT)];
        jj = [max(j(k)-1,1) j(k) min(j(k)+1,File(f).NumNT)];
        kk = find(E(i(k),:) == 1);
        N1 = File(f).NT(i(k));
        N2 = File(f).NT(j(k));
        zShowInteractionTable(File(f),[ii jj kk]);
        fprintf('Removing ncWW between %s %s%s and %s%s\n', File(f).Filename, N1.Base, N1.Number, N2.Base, N2.Number);
        figure(1)
        clf
        zDisplayNT(File(f),[ii jj kk]);
        drawnow
        pause
      end

    end
  end    

  if Verbose > 0 && length(i) > 0,
    fprintf('Removed %4d ncWW from %4s, which is %7.2f%% of the original number\n', r, File(f).Filename, 100*r/length(i));
  end

  % ----------------------------- remove multiple basepairs

  [i,j,e] = find(File(f).Edge);
  k = find(abs(e) < 2);             % cWW basepairs only
  i = i(k);
  j = j(k);
  e = e(k);

  u = zeros(size(e));

  for k = 1:length(i),
    t = zEdgeText(e(k));
    switch upper(t(2))
    case   'W', u(k) = 1;     % WC edge
    case   'H', u(k) = 2;     % Hoogsteen edge
    case   'S', u(k) = 3;     % sugar edge
    end
  end

  w = sparse(i,j,u);

  for b = 1:3,                          % which edge we are checking
    i = find(sum(w'==b) > 1);           % nucleotides using an edge twice
    m = [];                             % keep track of worst discrepancy
    for a = 1:length(i),                % loop through double-edge nucleotides
      j = find(w(i(a),:) == b);         % nucleotides that i(a) interacts with
      NT1 = File(f).NT(i(a));
      d = [];
      for c = 1:length(j),
        NT2 = File(f).NT(j(c));
        d(c) = zDistanceToExemplar(Exemplar,NT1,NT2,fix(File(f).Edge(i(a),j(c))));

        if b == 1 && NT1.Code == 2 && NT2.Code == 2,  % lousy CC pairs
          d(c) = 1;
        end

        if Verbose > 0,
          fprintf('Pair %s %s%5s_%s - %s%5s_%s %s %4.1f distance %7.4f to exemplar\n', File(f).Filename, NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File(f).Edge(i(a),j(c))), full(File(f).Edge(i(a),j(c))), d(c));
        end
      end
      m(a) = max(d);
    end

    [y,j] = sort(-m);                   % worst instances first
    i = i(j);                           % re-order double-edge nucleotides

    for a = 1:length(i),
      j = find(w(i(a),:) == b);         % nucleotides that i(a) interacts with
      if length(j) > 1,                 % we might have already handled this

        NT1 = File(f).NT(i(a));
        if Verbose > 2,
          clf
          VP.Sugar = 1;
%          zDisplayNT(File(f),[i(a) j i(a)-1 i(a)+1 min(j)-1 max(j) + 1],VP);
          zDisplayNT(File(f),[i(a) j],VP);
        end
        d = [];
        for c = 1:length(j),
          NT2 = File(f).NT(j(c));
          d(c) = zDistanceToExemplar(Exemplar,NT1,NT2,fix(File(f).Edge(i(a),j(c))));
                                        % distance to exemplar of main class

          if b == 1 && NT1.Code == 2 && NT2.Code == 2,  % lousy CC pairs
            d(c) = 1;
          end

        end
        for c = 1:length(j),
          if d(c) > min(d),
            if Verbose > 1,
              NT2 = File(f).NT(j(c));
              fprintf('Removing  %s%5s_%s - %s%5s_%s %s distance %7.4f to exemplar\n', NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File(f).Edge(i(a),j(c))), d(c));

            end
            if Verbose > 2,
              pause
            end

            w(i(a),j(c)) = 0;
            w(j(c),i(a)) = 0;

            e = File(f).Edge(i(a),j(c));
            if e > 0 && e < 100,
              e = e + 100;
            elseif e < 0 && e > -100,
              e = e - 100;
            end

            File(f).Edge(i(a),j(c)) = e;          % remove the worse basepair(s)
            File(f).Edge(j(c),i(a)) = -e;

            if Verbose > 1,
              fprintf('Now it is %s%5s_%s - %s%5s_%s %s distance %7.4f to exemplar\n', NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File(f).Edge(i(a),j(c))), d(c));
            end

          end
        end
        if Verbose > 0,
%          fprintf('\n');
        end
      end
    end
  end  
end

if Verbose > 0,
  fprintf('Removed %4d ncWW, which is %7.2f%% of the original number\n', rr, 100*rr/nn);
end
