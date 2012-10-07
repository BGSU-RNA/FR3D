% zInteractionRange(File) returns a sparse matrix S which indicates
% how far a pairwise interaction is from being a local interaction

function [File] = zInteractionRange(File, Verbose)

if nargin < 2,
  Verbose = 1;
end

for f = 1:length(File),

if length(File.NT) > 1,

C = fix(abs(File(f).Edge))+0.001*(fix(abs(File(f).BasePhosphate+File(f).BasePhosphate')));                  % matrix of interactions, including base-phosphate
E = fix(abs(File(f).Edge)); % matrix of base-base interactions

for i = 1:length(C(:,1)),
  C(i,i) = 0;               % remove phosphate self interactions
end

C = triu(C);                % sparse matrix of pairwise interactions

[i,j,c] = find(C);          % locations of all pairwise interactions

c = fix(c);                 % remove decimal from base-phosphate

k = find((j > i));          % ignore self interactions

i = i(k);
j = j(k);
c = c(k);

[y,k] = sort(j-i);          % order by distance from diagonal

i = i(k);
j = j(k);
c = c(k);

d = zeros(size(c));                              % current interaction range

for m = 1:length(i),                             % go through all pairs once
  if (c(m) == 1) && (d(m) == 0),                 % cWW pair, not pseudoknotted

    p =     (i <= i(m)) .* (j >= i(m)) .* (j <= j(m));
    p = p + (i >= i(m)) .* (i <= j(m)) .* (j >= j(m));
    q = find(p);                      % indices of pairs that overlap i(m),j(m)

%    d(q) = max(d(q), abs(i(m)-i(q))+abs(j(m)-j(q))); 
                            % increase interaction range for these pairs
                            % to account for how far they cross i(m),j(m)

    d(q) = max([d(q) min([abs(i(m)-i(q)) abs(j(m)-j(q))],[],2)],[],2); 
                            % for all nested cWW pairs crossed, keep track
                            % of how many nucleotides this pair would need to
                            % move over to get uncrossed, keep the maximum

%    d(q) = d(q) + 1;        % count number of nested cWW's crossed

    d(m) = 0;               % don't inadvertently make this non-nested!

    if Verbose > 2,
      figure(2)
      clf
      plot(j(q),i(q),'.');
      hold on
      plot(j(m),i(m),'r.');
      axis([1 max(i) 1 max(i)]);
      axis ij
      drawnow
      if Verbose > 3,
        pause
      end
    end

  end
end

d = min(d,abs(i-j));                % adjust for local in-strand interactions

i = [i; length(File(f).NT)];
j = [j; length(File(f).NT)];
d = [d; 0];
c = [c; 0];

S = sparse(i,j,d);
S = S + S';
File(f).Range = S;

if Verbose > 1,
  fprintf('%s has %d basepairs, of which %d are local.\n', File.Filename, full(sum(sum((C > 0) .* (C < 15)))), full(sum(sum(((S<=10).*C > 0) .* ((S<=10).*C < 15)))));
  for m = 1:14,
    all = length(find(c==m));
    local = length(find((c==m).*(d<=10)));
    longrange = all - local;
    fprintf('%4s %3d total, %3d local, %3d long-range\n', zEdgeText(m), all, local, longrange);
  end

  fprintf('Long-range cWW interactions:\n');
  q = find((c==1).*(d>10));
  [y,h] = sort(i(q));
  q = q(h);
  for m = 1:length(q),
    fprintf('%s%s with %s%s\n', File.NT(i(q(m))).Base,File.NT(i(q(m))).Number, File.NT(j(q(m))).Base, File.NT(j(q(m))).Number);
  end
end

else

  File.Range = [];

end

end
