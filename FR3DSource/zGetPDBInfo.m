% zGetPDBInfo(File) reads a data file about PDB files and extracts
% information about the current molecule stored in File

% It uses PDBInfo.mat, which needs to be updated whenver new PDB data is
% downloaded from NDB/PDB

function [File] = zGetPDBInfo(File,n,t)

if nargin < 2,
  load(['FR3DSource' filesep 'PDBInfo.mat'],'n','t','-mat');
end

PDBNames = lower(t(:,1));              % convert PDB filenames to lowercase

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

