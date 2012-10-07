
function [File] = zBackboneContinuity(File)

for f = 1:length(File),

  File(f).Covalent = sparse(zeros(File(f).NumNT,File(f).NumNT));

  d = [];
  g = [];

  if isempty(File(f).Distance),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angstroms
  end

  for i = 2:length(File(f).NT),
    a = File(f).NT(i).Sugar(10,:);      % phosphorus of current NT
    b = File(f).NT(i-1).Sugar(5,:);     % O3 of previous NT
    c = File(f).NT(i-1).Sugar(3,:);     % O2 of previous NT

    x = norm(a-b);
    y = norm(a-c);

    if x < 2  && x < y,
      File(f).Covalent(i-1,i) = 1;
      File(f).Covalent(i,i-1) = -1;
    elseif y < 2,
      File(f).Covalent(i-1,i) = 2;
      File(f).Covalent(i,i-1) = -2;
    else 
      j = find(File(f).Distance(i-1,:));     % nucleotides somewhat near i-1
      for k = 1:length(j),
        a = File(f).NT(j(k)).Sugar(10,:);      % phosphorus of current NT
        x = norm(a-b);
        y = norm(a-c);

        if x < 2  && x < y,
          File(f).Covalent(i-1,j(k)) = 1;
          File(f).Covalent(j(k),i-1) = -1;
fprintf('Covalent connection c35 between nucleotides not in order in the PDB file.  %s %s %s\n', File(f).Filename, File(f).NT(i-1).Number, File(f).NT(j(k)).Number);
        elseif y < 2,
          File(f).Covalent(i-1,j(k)) = 2;
          File(f).Covalent(j(k),i-1) = -2;

fprintf('Covalent connection c25 between nucleotides not in order in the PDB file.  %s %s %s\n', File(f).Filename, File(f).NT(i-1).Number, File(f).NT(j(k)).Number);
        end
      end
    end

%    d = [d; norm(a-b)];
%    g = [g; norm(a-c)];


  end
end

return

min(d)
max(d)

e = sort(d);

e(end-10:end)

figure(2)
clf
i = find(d < 2);
hist(d(i),30)

figure(3)
clf
i = find(d < 2);
plot(d(i),g(i),'.')

figure(1)
clf
VP.Sugar = 1;
VP.LabelSugar = 10;
VP.LabelBase = 10;
zDisplayNT(File(1),[10 11],VP);
axis square
