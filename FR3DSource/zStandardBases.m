% zStandardBases reads optimized bases, then centers and rotates them

% Variables returned to the rest of the programs:
% 
% StandardLoc      16 x 3 x 4 matrix of standard atom locations, all 4 bases
% AtomNames     16 x 4 cell array of names of atoms in standard bases
% Lim               2 x 4 matrix telling how many base atoms, total atoms
%                         (including hydrogen) the 4 bases have
%
% A_Stand, C_Stand, G_Stand, U_Stand: rotated versions of standard bases

% ------------------------------------------------------------------------
[A_Atoms, zX, zY, zZ] = zStandardBase('A');

Q = [zX zY zZ];
Q = Q - ones(length(zX),1)*Q(7,:);    % center base atoms at N9
v = Q(13,:) - Q(1,:);
c = v(2)/norm(v);
s = v(1)/norm(v);
Q = -Q * [[c s 0]; [-s c 0]; [0 0 1]]; % rotate 
A_Stand = Q;

% ------------------------------------------------------------------------
[C_Atoms, zX, zY, zZ] = zStandardBase('C');

Q = [zX -zY zZ];
Q = Q - ones(length(zX),1)*Q(2,:);
v = Q(9,:) - Q(1,:);
c = v(2)/norm(v);
s = v(1)/norm(v);
Q = -Q * [[c s 0]; [-s c 0]; [0 0 1]];
C_Stand = Q;

% ------------------------------------------------------------------------
[G_Atoms, zX, zY, zZ] = zStandardBase('G');

Q = [zX -zY zZ];
Q = Q - ones(length(zX),1)*Q(10,:);
v = Q(14,:) - Q(1,:);
c = v(2)/norm(v);
s = v(1)/norm(v);
Q = -Q * [[c s 0]; [-s c 0]; [0 0 1]];
G_Stand = Q;

% ------------------------------------------------------------------------
[U_Atoms, zX, zY, zZ] = zStandardBase('U');

Q = [zX -zY zZ];
Q = Q - ones(length(zX),1)*Q(3,:);
v = Q(10,:) - Q(1,:);
c = v(2)/norm(v);
s = v(1)/norm(v);
Q = -Q * [[c s 0]; [-s c 0]; [0 0 1]];
U_Stand = Q;

%----------------------------------------------------------------------------
StandardLoc = zeros(16,3,4);

% Center bases at glycosidic atom

StandardLoc(1:15,:,1)= A_Stand-ones(size(A_Stand(:,1)))*A_Stand(1,:);
StandardLoc(1:13,:,2)= C_Stand-ones(size(C_Stand(:,1)))*C_Stand(1,:);
StandardLoc(1:16,:,3)= G_Stand-ones(size(G_Stand(:,1)))*G_Stand(1,:);
StandardLoc(1:12,:,4)= U_Stand-ones(size(U_Stand(:,1)))*U_Stand(1,:);

Lim(1,:) = [10 8 11 8];       % number of base atoms, excluding hydrogen
Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen
Lim(3,:) = [13  9 14 10];     % locations of fictitious hydrogens

AtomNames = cell(16,4);

AtomNames(1:Lim(2,1),1) = A_Atoms;   % provide a list of atom names
AtomNames(1:Lim(2,2),2) = C_Atoms;
AtomNames(1:Lim(2,3),3) = G_Atoms;
AtomNames(1:Lim(2,4),4) = U_Atoms;


% code to output standard base locations as Python code

if 0 > 1,
	CenterLoc(1:15,:,1)= A_Stand-ones(15,1)*mean(A_Stand(1:10,:));
	CenterLoc(1:13,:,2)= C_Stand-ones(13,1)*mean(C_Stand(1:8,:));
	CenterLoc(1:16,:,3)= G_Stand-ones(16,1)*mean(G_Stand(1:11,:));
	CenterLoc(1:12,:,4)= U_Stand-ones(12,1)*mean(U_Stand(1:8,:));
  for b = 1:4,
  	sum(CenterLoc(1:Lim(1,b),:,b))
  end
	RNA = 'ACGU';
	fprintf('# Quantum mechanics optimized base atom locations, by Jiri Sponer\n');
  for b = 1:4,
  	for a = 1:Lim(2,b),
  	  fprintf('RNAbasecoordinates[''%s''][% 5s] = [% 10.6f, % 10.6f, % 10.6f]\n', RNA(b), ['''' AtomNames{a,b} ''''], CenterLoc(a,1,b), CenterLoc(a,2,b), CenterLoc(a,3,b));
  	end
  end
end