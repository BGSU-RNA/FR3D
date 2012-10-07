
function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, NTLETTER, CHAIN, NTNUMBER, P,OCC,TEMP] = zReadPDBTextRead(Filename)

Readable = 1;

try

  [A, B, C, E, F, G, X, Y, Z, OCC, TEMP] ...
   = textread(strcat(Filename,'.pdb'),'%6s%5d  %4c%4s %1s%5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');

catch

  try
    [A, B, C, E, G, X, Y, Z, OCC, TEMP] ...
     = textread(strcat(Filename,'.pdb'),'%6s%5d  %4c%4s %5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');
    for i=1:length(A),
      F{i} = '1';                % invent a chain number
    end
  catch

    fprintf('Unable to read the PDB file %s\n',Filename);
    Readable = 0;

  end
end

if Readable > 0,

[s,t] = size(C);

NoChain = 0;

for i=1:s,
  if ~isempty(strfind(G{i},'.')),
    NoChain = 1;            % when no chain info, NTNUMBER is read from X column
  end
end

CHAIN = cell(s,1);

if NoChain == 1,
  [A, B, C, E, G, X, Y, Z, OCC, TEMP] ...
   = textread(strcat(Filename,'.pdb'),'%6s%5d  %4c%4s %5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');
  for i=1:s,
    CHAIN{i,1} = '1';
  end
else
  CHAIN      = [F];
end

ATOM_TYPE  = [A];
ATOMNUMBER = [B];
NTLETTER   = [E];
NTNUMBER   = [G];
P          = [X Y Z];

ATOMNAME = cell(s,1);
VERSION  = cell(s,1);

for i=1:s,
  ATOMNAME{i,1} = deblank(C(i,1:3));
  VERSION{i,1}  = C(i,4);
end

%fprintf('Read  %s\n', [Filename '_Atoms.pdb']);

else

ATOM_TYPE = [];
ATOMNUMBER = [];
NTLETTER = [];
NTNUMBER = [];
P        = [];
ATOMNAME = [];
VERSION  = [];
OCC      = [];
TEMP     = [];

end
