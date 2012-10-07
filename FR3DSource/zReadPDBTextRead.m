% some PDB files have some atoms with two "versions" A and B.

function [ATOM_TYPE, ATOMNUMBER, ATOMNAME, VERSION, NTLETTER, CHAIN, NTNUMBER, P,OCC,TEMP,Readable] = zReadPDBTextRead(Filename)

Readable = 1;

try

disp('zReadPDBTextRead method 1');

  [A, B, C, E, F, G, X, Y, Z, OCC, TEMP] = textread(Filename,'%5s%6d  %4c%3c %1c%4s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');

%  [A, B, C, E, F, G, X, Y, Z, OCC, TEMP] = textread(Filename,'%4s%7d%5s%4s%2s%5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');

catch

  try


disp('zReadPDBTextRead method 2');

    [A, B, C, E, G, X, Y, Z, OCC, TEMP] ...
     = textread(Filename,'%5s%6d  %4c%4s %5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');

    for i=1:length(A),
      F{i} = 'Z';                % invent a chain number
    end
  catch

    try

disp('zReadPDBTextRead method 3');


      [A, B, C, F, G, X, Y, Z, OCC, TEMP] = textread(Filename,'%4c%7d%8c%2s%5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');

C

E = C(:,(end-3):end);

E

    catch

      Readable = 0;

    end
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
   = textread(Filename,'%5s%6d  %4c%4s %5s  %8.3f%8.3f%8.3f%6.2f%6.2f%*[^\n]');
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
CHAIN    = [];

end
