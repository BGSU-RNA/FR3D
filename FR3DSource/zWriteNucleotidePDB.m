% zWriteNucleotidePDB(fid,NT,a,c,R,Sh) writes the atoms of nucleotide NT, with
% hydrogens, to the current file id.
% If specified, it rotates by R, shifts by sh, puts each c set off from
% the others by 20*c Angstroms

function [a] = zWriteNucleotidePDB(fid,NT,a,c,R,Sh)

if nargin < 4,
  c = 0;
end

if nargin < 5,
  c = 0;
  R = eye(3);
  Sh = [0 0 0];
end

x = mod(c,30);                       % shift along x axis
y = mod(fix(c/30),30);               % shift along y axis
z = mod(fix(c/900),900);             % shift along z axis

Lim(1,:) = [10 8 11 8 Inf Inf];     % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12 Inf Inf];   % total number of atoms, including hydrogen
Lim(3,:) = [13  9 14 10 Inf Inf];   % locations of fictitious hydrogens

A = {' N9' ' C4' ' N3' ' N1' ' C6' ' N6' ' C8' ' C5' ' C2' ' N7' ' H2' ' H8' ' H9' '1H6' '2H6'};
C = {' N1' ' C2' ' O2' ' N3' ' C4' ' N4' ' C6' ' C5' ' H1' ' H6' ' H5' '1H4' '2H4'};
G = {' N9' ' C4' ' N3' ' N1' ' C6' ' O6' ' C8' ' C5' ' C2' ' N7' ' N2' ' H1' ' H8' ' H9' '1H2' '2H2'};
U = {' N1' ' C2' ' O2' ' N3' ' C4' ' O4' ' C6' ' C5' ' H5' ' H1' ' H3' ' H6'};
AA= {'N' 'C' 'N' 'N' 'C' 'N' 'C' 'C' 'C' 'N' 'H' 'H' 'H' 'H' 'H'};
CA= {'N' 'C' 'O' 'N' 'C' 'N' 'C' 'C' 'H' 'H' 'H' 'H' 'H'};
GA= {'N' 'C' 'N' 'N' 'C' 'O' 'C' 'C' 'C' 'N' 'N' 'H' 'H' 'H' 'H' 'H'};
UA= {'N' 'C' 'O' 'N' 'C' 'O' 'C' 'C' 'H' 'H' 'H' 'H'};
S = {' C1*' ' C2*' ' O2*' ' C3*' ' O3*' ' C4*' ' O4*' ' C5*' ' O5*' ' P  ' ' OP1' ' OP2'};
SA= {'C' 'C' 'O' 'C' 'O' 'C' 'O' 'C' 'O' 'P' 'O' 'O'};

A = strrep(A,'*','''');
C = strrep(C,'*','''');
G = strrep(G,'*','''');
U = strrep(U,'*','''');
S = strrep(S,'*','''');

%      1     2     3     4     5     6     7     8     9     10  11    12
SugarReorder = [10 11 12 9 8 6 7 4 5 2 3 1];

  for k = 1:min(12,length(NT.Sugar(:,1))),         % loop through sugar atoms
    j = SugarReorder(k);
    fprintf(fid, 'ATOM %6d', a);
    fprintf(fid,' %3s', S{j});
    fprintf(fid, ' %3s', NT.Base);
    fprintf(fid, ' %1s', NT.Chain);
    fprintf(fid, '%4s    ', NT.Number);
    L = (NT.Sugar(j,:) - Sh)*R;
    L = L + 30*[x y z];
    fprintf(fid, '%8.3f', L);
    fprintf(fid, '%6.2f', 1);
    fprintf(fid, '%6.2f', 99.99);
    fprintf(fid, '           %c\n',SA{j});
    a = a + 1;
  end
  s = length(NT.Fit(:,1));

  for j = 1:min(Lim(2,NT.Code),s),             % loop through base atoms and H
    if j ~= Lim(3,NT.Code) && NT.Code ~= 6,
      fprintf(fid, 'ATOM %6d', a);
      switch NT.Code,
        case 1, fprintf(fid,' %3s', A{j});
        case 2, fprintf(fid,' %3s', C{j});
        case 3, fprintf(fid,' %3s', G{j});
        case 4, fprintf(fid,' %3s', U{j});
        case 5, fprintf(fid,' %3s', NT.AtomName{j});  % modified
        case 6, fprintf(fid,' %3s', NT.AtomName{j});  % missing base atom(s)
      end
      fprintf(fid, '  %3s', NT.Base);
      fprintf(fid, ' %1s',  NT.Chain);
      fprintf(fid, '%4s    ',   NT.Number);
      L = (NT.Fit(j,:) - Sh)*R;
      L = L + 30*[x y z];                        % shift 20 Angstroms
      fprintf(fid, '%8.3f', L);                  % write atom location
      fprintf(fid, '%6.2f', 1);
      fprintf(fid, '%6.2f', 99.99);
      switch NT.Code,
        case 1, fprintf(fid,'           %c\n', AA{j});
        case 2, fprintf(fid,'           %c\n', CA{j});
        case 3, fprintf(fid,'           %c\n', GA{j});
        case 4, fprintf(fid,'           %c\n', UA{j});
        case 5, fprintf(fid,'           %c\n', NT.AtomName{j});  % modified
        case 6, fprintf(fid,'           %c\n', NT.AtomName{j});  % missing base atom(s)
      end
      a = a + 1;
    end
  end
