% zGetPDBInfo(File) attempts to retrieve information from PDB if it is not already present

function [File] = zGetPDBInfo(File)

t = '';

PDBID = File.PDBID;
PDBID = strrep(PDBID,'-CIF','');              % strip off this identifier of the source of the file
PDBID = strrep(PDBID,'.mat','');

attempts = 0;
while attempts < 30,
  attempts = attempts + 1;
  try
    t = urlread(['http://www.pdb.org/pdb/rest/customReport.xml?pdbids=' PDBID '&customReportColumns=source,resolution,experimentalTechnique,structureTitle,releaseDate']);
  catch
    pause(1)
    fprintf('zGetPDBInfo:  Having trouble reading %s data from PDB; attempt %d\n',PDBID,attempts);
    t
  end
end

i = strfind(t,'<dimStructure.resolution>');
j = strfind(t,'</dimStructure.resolution>')-1;
r = strrep(t(i(1):j(1)),'<dimStructure.resolution>','');
File.Info.Resolution = str2num(r);
if isempty(File.Info.Resolution),
  File.Info.Resolution = NaN;
end

i = strfind(t,'<dimStructure.experimentalTechnique>');
j = strfind(t,'</dimStructure.experimentalTechnique>')-1;
r = strrep(t(i(1):j(1)),'<dimStructure.experimentalTechnique>','');
File.Info.ExpTechnique = r;

i = strfind(t,'<dimStructure.structureTitle>');
j = strfind(t,'</dimStructure.structureTitle>')-1;
r = strrep(t(i(1):j(1)),'<dimStructure.structureTitle>','');
File.Info.Descriptor = r;

i = strfind(t,'<dimStructure.releaseDate>');
j = strfind(t,'</dimStructure.releaseDate>')-1;
r = strrep(t(i(1):j(1)),'<dimStructure.releaseDate>','');
File.Info.ReleaseDate = r;

i = strfind(t,'<dimStructure.releaseDate>');
j = strfind(t,'</dimStructure.releaseDate>')-1;
r = strrep(t(i(1):j(1)),'<dimStructure.releaseDate>','');
File.Info.Author = r;

i = strfind(t,'<dimEntity.chainId>') + 19;
j = strfind(t,'</dimEntity.chainId>') - 1;
m = strfind(t,'<dimEntity.source>') + 18;
n = strfind(t,'</dimEntity.source>') - 1;

if isfield(File,'LongestChain') && ~isempty(File.LongestChain) && ~isempty(File.NT),

  LC = File.NT(File.LongestChain(1)).Chain;          % identifier of longest chain

  for k = 1:length(i),
    if strcmp(LC,t(i(k):j(k))),
      File.Info.Source = t(m(k):n(k));
    end
  end
else
  File.Info.Source = t(i(1):j(1));
end

return

% i = strmatch(lower(File.Filename),PDBNames);

i = find(ismember(PDBNames, lower(File.Filename)));

if ~isempty(i),
  File.Info.Resolution  = n(i(1),1);
  File.Info.Descriptor  = t{i(1),2};
  File.Info.ExpTechnique= t{i(1),3};
  File.Info.ReleaseDate = t{i(1),4};
  File.Info.Author      = t{i(1),5};
  File.Info.Keywords    = t{i(1),6};
  File.Info.Source      = t{i(1),8};

  if ((n(i(1),2) ~= length(File.NT)) || isempty(t{i(1),9}) || n(i(1),3) == 0) && (length(File.NT) > 0) ,
    E  = abs(triu(File.Edge));
    P = zSparseRange(E,0,16);
    n(i(1),3) = nnz(P); % number of pairs
    n(i(1),2) = length(File.NT);
    t{i(1),9} = cat(2,File.NT.Base);
%    save(['FR3DSource' filesep 'PDBInfo.mat'],'n','t'); % Matlab version 7
%    don't push data back to PDBInfo.mat; rely on zUpdatePDBDatabase to do that
  end

else
  File.Info.Resolution  = [];
  File.Info.Descriptor  = '';
  File.Info.ExpTechnique= '';
  File.Info.ReleaseDate = '';
  File.Info.Author      = '';
  File.Info.Keywords    = '';
  File.Info.Source      = '';
end

