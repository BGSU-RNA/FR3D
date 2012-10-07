% zInteractionRange(File) returns a sparse matrix S which indicates
% how far a pairwise interaction is from being a local interaction

function [File] = zInteractionRange(File, Verbose)

for f = 1:length(File),

if length(File.NT) > 1,

C = fix(abs(File(f).Edge))+0.01*(fix(abs(File(f).BasePhosphate+File(f).BasePhosphate')));
E = fix(abs(File(f).Edge));

for i = 1:length(C(:,1)),
  C(i,i) = 0;               % remove phosphate self interactions
end

if Verbose > 3,
  figure(1)
  clf
  zNussinovPlot(File, triu(E.*(abs(E)<20)), E.*(E<20));
  title('Locations of basepair interactions of all types');
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

d = zeros(size(c));                              % current distance from local

figure(2)
spy(triu(E.*(abs(E)<20)));

for m = 1:length(i),
  if (c(m) == 1) && (d(m) <= 10),                % pair, not pseudoknot
    p =     (i <= i(m)) .* (j >= i(m)) .* (j <= j(m)); % pairs whose interaction range are affected by the current basepair
    p = p + (i >= i(m)) .* (i <= j(m)) .* (j >= j(m)); % more such pairs
    q = find(p);                                 % indices of these pairs

    if c(m) == 1,                      % current basepair is cWW
      a = 1;                           % record Manhattan distance to i(m),j(m)
    else
      a = 0.5;                         % use half the distance
    end

    d(q) = max(d(q), a*(abs(i(m)-i(q))+abs(j(m)-j(q)))); % Manhattan distance

    if Verbose > 2,
      figure(3)
      clf
      plot(j(q),i(q),'.');             
      hold on
      plot(j(m),i(m),'r.');            % current basepair
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

if Verbose > 1,
  figure(3)
  clf
  S = sparse(i,j,d);
  zNussinovPlot(File, triu((E == 1)), 1+2*(S > 0));
  title(['cis Watson-Crick interactions in ' File(f).Filename]);

  figure(4)
  clf
  B = E .* (E > 0) .* (E < 24);                 % pairs and stacks
  zNussinovPlot(File,triu(B), (B==1).*(S==0) + 2*(B>1).*(B<14).*(S==0) + 3*(B==1).*(S>0) + 4*(B > 1).*(B < 14) .*(S>0) + 5*(B > 20) .* (B < 25) .* (S > 10),0.8);
%  title(['All basepairs in ' File(f).Filename]);

end

S = sparse(i,j,d);
S = S + S';
File(f).Range = S;

if Verbose > 0,
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
