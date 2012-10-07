% zBaseRiboseTable(File,NTList) displays a table of interactions
% among the nucleotides NTList.  If passed, it displays the number Disc

function [void] = zBaseRiboseTable(File,NTList,Disc)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

%Indices = sort(Indices);

fprintf('  File %s',File.Filename);         % display filename
fprintf(' Chain ');
for j=1:length(Indices),
  fprintf('%s',File.NT(Indices(j)).Chain);
end
fprintf('\n');

if nargin == 3,
  fprintf('%6.4f',Disc);                 % display discrepancy if passed
else
  fprintf('       ');
end

for j=1:length(Indices),
  fprintf('%6s',[File.NT(Indices(j)).Base File.NT(Indices(j)).Number]);
end
fprintf('\n');

Config = {'(A)' , '(S)'};

for i=1:length(Indices),
  fprintf('%6s',[File.NT(Indices(i)).Base File.NT(Indices(i)).Number]);
  for j=1:length(Indices),
    fprintf('%6s', zBaseRiboseText(File.BaseRibose(Indices(i),Indices(j))));
  end
  fprintf('\n');
end

drawnow
