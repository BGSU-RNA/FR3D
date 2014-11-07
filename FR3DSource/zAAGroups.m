% zAAGroups(Unit,Atom) returns two lists of indices into the cell array Atom, depending on the Unit, which is a 3-letter identifier for amino acids.
% GroupOne is the "business end" of the amino acid sidechain and GroupTwo typically contains the backbone and sometimes additional atoms

function [GroupOne,GroupTwo,GroupThree] = zAAGroups(Unit,Atom)

GroupTwoAtoms = {};

switch Unit
case 'ARG'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};           % backbone atoms to be colored blue
  GroupTwoAtoms = {'CB','CG','CD'};                      % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'NE','CZ','NH1','CZ','NH2'};    % side chain or business end

case 'LYS'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};           % backbone atoms to be colored blue
  GroupTwoAtoms = {'CB','CG','CD','CE'};               % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'NZ'};                                 % side chain or business end

case 'HIS'
   GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
   GroupTwoAtoms = {'CB'};                                                           % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD2','NE2','CE1','ND1','CG'};                  % side chain or business end

case 'GLN'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupTwoAtoms = {'CB'};                                                          % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD','OE1','CD','NE2'};                           % side chain or business end

case 'ASN'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupThreeAtoms = {'CB','CG','OD1','CG','ND2'};                           % side chain or business end

case 'PHE'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupTwoAtoms = {'CB'};                                                          % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD1','CE1','CZ','CE2','CD2','CG'};        % side chain or business end

case 'TYR'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                                      % backbone atoms to be colored blue 
  GroupTwoAtoms = {'CB'};                                                                              % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD1','CE1','CZ','OH','CZ','CE2','CD2','CG'};               % side chain or business end

case 'TRP'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                                              % backbone atoms to be colored blue 
  GroupTwoAtoms = {'CB'};                                                                                      % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD1','NE1','CE2','CD2','CE3','CZ3','CH2','CZ2','CE2'};          % side chain or business end

case 'ASP'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                    % backbone atoms to be colored blue 
  GroupThreeAtoms = {'CB','CG','OD1','CG','OD2'};                           % side chain or business end

case 'GLU'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupTwoAtoms = {'CB'};                                                          % list of atoms to join with lines, in this order
  GroupThreeAtoms = {'CG','CD','OE1','CD','OE2'};                          % side chain or business end

case 'VAL'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupThreeAtoms = {'CB','CG1','CB','CG2'};                                 % side chain or business end

case 'ALA'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue 
  GroupThreeAtoms = {'CB'};                                 % side chain or business end

case 'LEU'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue   
  GroupThreeAtoms = {'CB','CG','CD1','CG','CD2'};                           % side chain or business end

case 'ILE'
  GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue    
  GroupThreeAtoms = {'CB','CG1','CB','CG2','CD1'};                         % side chain or business end

case 'MET'
 GroupOneAtoms = {'N','CA','C','O','C','CA'};                                   % backbone atoms to be colored blue      
 GroupThreeAtoms = {'CB','CG','SD','CE'};          % list of atoms to join with lines, in this order
 
case 'GLY'
  GroupOneAtoms = {'N','CA','C','O'};   % backbone atoms to be colored blue 
  GroupThreeAtoms = {'CB'};             % but sometimes GLY doesn't have this atom?

otherwise
  GroupOneAtoms = {'N','CA','C','O'};
  GroupThreeAtoms = setdiff(Atom,GroupTwoAtoms);
  fprintf('Unknown amino acid %s\n',Unit);
end

GroupOne = [];
for j = 1:length(GroupOneAtoms),
  a = find(ismember(Atom,GroupOneAtoms{j}));
  if ~isempty(a),
    GroupOne = [GroupOne a(1)];
  end
end

GroupTwo = [];
for j = 1:length(GroupTwoAtoms),
  a = find(ismember(Atom,GroupTwoAtoms{j}));
  if ~isempty(a),
    GroupTwo = [GroupTwo a(1)];
  end
end

GroupThree = [];
for j = 1:length(GroupThreeAtoms),
  a = find(ismember(Atom,GroupThreeAtoms{j}));
  if ~isempty(a),
    GroupThree = [GroupThree a(1)];
  end
end

Extras = setdiff(Atom,union(GroupOneAtoms,union(GroupTwoAtoms,GroupThreeAtoms)));
if ~isempty(Extras),
  fprintf('Unexpected atoms found in amino acid %s\n',Unit);
  Extras
end

if 0 > 1,
  Unit
  Atom
  GroupOneAtoms
  GroupOne
  GroupTwoAtoms
  GroupTwo

  pause
end
