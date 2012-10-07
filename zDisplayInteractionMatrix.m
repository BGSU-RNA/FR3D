% zDisplayInteractionMatrix(File,NTList) displays a simple graph of
% the matrix of interactions among bases in NTList
% File and NTList may be specified as in zDisplayNT

function [void] = zDisplayInteractionMatrix(File,NTList)

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

figure(1)
clf

spy((abs(File.Inter(Indices,Indices)) > 0) .* (File.Inter(Indices,Indices)<30));
