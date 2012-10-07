


% check center to center distance between actual pairs

d = [];
c = 1;

for f = 1:length(File),
  [i,j] = find((abs(File(f).Edge) > 0) .* (abs(File(f).Edge) < 30));
  for k = 1:length(i),
    d = [d; File(f).Distance(i(k),j(k))];
  end
end

hist(d,30);