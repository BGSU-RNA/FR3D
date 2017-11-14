% zLoadIFE(Filename,ReadCode,Verbose) parses Filename to look for chain or model or multiple IFEs and loads just those
% File = zLoadIFE('4V4Q|1|CA')
% File = zLoadIFE('4V4Q|1|CA+4V4Q|1|DB')

function [File] = zLoadIFE(Filename,ReadCode,Verbose,File)

if nargin < 2,
  ReadCode = 0;
end

if nargin < 3,
  Verbose = 0;
end

IFE = zStringSplit(Filename,'+');             % split multiple chains in the IFE
a = zStringSplit(IFE{1},'|');                 % split by divider for filename, model, etc.

if nargin < 4,
  File = zGetNTData(a{1},ReadCode,Verbose);     % load the whole 3D structure file
end

if length(a) > 1,                              % a model or chain is specified
  keep = zeros(1,length(File.NT));             % 0 means not to keep a nucleotide
  for k = 1:length(IFE),
    for n = 1:length(File.NT),
      if ~isempty(strfind(File.NT(n).ID,IFE{k})),  % check each IFE part against each nucleotide ID
        keep(n) = 1;
      end
    end
  end

  indices = find(keep);

  if Verbose > 0,
    numnt = length(indices);
    fprintf('Keeping %5d nucleotides from %s\n',length(indices),IFE{1});
  end

  try
    File = zSubFile(File,indices);
  catch
    fprintf('zLoadIFE:  Something went wrong with making a subfile\n');
  end

%  File.IFE = IFE{1};            % save the IFE being saved here

end
