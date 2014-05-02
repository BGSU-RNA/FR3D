% zGetPDBInfo(File) attempts to retrieve information from PDB if it is not already present

function [File] = zGetPDBInfo(File)

t = urlread(['http://www.pdb.org/pdb/rest/customReport.xml?pdbids=' File.Filename '&customReportColumns=resolution,experimentalTechnique,structureTitle,releaseDate']);

i = strfind(t,'<dimStructure.resolution>');
j = strfind(t,'</dimStructure.resolution>')-1;
r = strrep(t(i:j),'<dimStructure.resolution>','');
File.Info.Resolution = str2num(r);

i = strfind(t,'<dimStructure.experimentalTechnique>');
j = strfind(t,'</dimStructure.experimentalTechnique>')-1;
r = strrep(t(i:j),'<dimStructure.experimentalTechnique>','');
File.Info.ExpTechnique = r;

i = strfind(t,'<dimStructure.structureTitle>');
j = strfind(t,'</dimStructure.structureTitle>')-1;
r = strrep(t(i:j),'<dimStructure.structureTitle>','');
File.Info.Descriptor = r;

i = strfind(t,'<dimStructure.releaseDate>');
j = strfind(t,'</dimStructure.releaseDate>')-1;
r = strrep(t(i:j),'<dimStructure.releaseDate>','');
File.Info.ReleaseDate = r;

File.Info

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

