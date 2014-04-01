% oExploreNT(File,NTList) is a general-purpose nucleotide plotting program.  It can be called in several ways, for example,
% oExploreNT('1s72',{'27','28','29'}), and it will load the
% named datafile and plot the nucleotides by nucleotide number.  Or, if
% data files have already been loaded, one can use oExploreNT(File(1),[32
% 34 35]) to plot the nucleotides in File(1) having indices 32,
% 34, and 35.
% One can also use ranges of nucleotide numbers, as in
% oExploreNT('1S72',{'2548:2555','2557','2559:2566'});

function [File] = oExploreNT(File,NTList)

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

if nargin == 1 || isempty(NTList),
  NTList = 1:length(File.NT);                  % display them all
end

% if NTList is a cell array of numbers, look up the indices

if strcmp(class(NTList),'char'),
  NTList = strrep(NTList,'/',' ');         % some people use this in papers
  NTList = {NTList};
end

if strcmp(class(NTList),'cell'),
  Indices = zIndexLookup(File,NTList);
else
  Indices = NTList;
end

TextNTList = [];
ChainList = [];

for j = 1:min(length(Indices),10),
  i = Indices(j);
  TextNTList = [TextNTList '_' File.NT(i).Base File.NT(i).Number];
  ChainList  = [ChainList File.NT(i).Chain];
end
TextNTList = [TextNTList '_' ChainList];

Search.File = File;
Search.Candidates = [Indices 1];
Search.SaveName = ['Explore_' File.Filename TextNTList];
Search.oExploreNT = 1;

xDisplayCandidates(File,Search);
