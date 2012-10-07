% xSpecifyQuery returns the description of a model motif
%
% The variable Query has several fields, most of which are optional:
%   Query.Description    a useful string, can be long
%   Query.Name           a short string, which will become part of a filename
%
% For all searches, these are required:
%   Query.SearchFiles    a cell array of PDB filenames or PDB lists
%
% For geometric searches, these are required:
%   Query.Filename       a string like 1s72, where the query motif is found
%   Query.NTList         a list of nucleotide numbers
%   Query.ChainList      a list of chain specifications, use if needed
%
% For non-geometric or mixed searches, these are optional:
%   Query.Edges          list of required basepairing or stacking interactions
%                        also specify allowed or disallowed pairs here
%   Query.Mask           a mask for which nucleotides to allow (defaults N)
%   Query.AngleWeight    weights to put on the angles (defaults 1)
%   Query.DistanceWeight weights to put on nucleotide distances (defaults 1)
%   Query.DiscCutoff     discrepancy cutoff D_0 (default 0.4)
%   Query.RelCutoff      relaxed cutoff D_1 (default Query.DiscCutoff)
%   Query.MaxDiff        maximum difference between sorted nucleotide indices
%   Query.MinDiff        minimum difference between sorted nucleotide indices
%
%   Query.Geometric      set to 0 to ignore geometry, only use screens. 
%                        Default is 1.
%   Query.ExcludeOverlap set to 1 to eliminate highly redundant motifs with
%     larger discrepancies; often, you get the same candidate with one or two
%     nucleotides different but much higher discrepancy.  Default is 1 when
%     Query.NumNT > 6.

function [Query] = xSpecifyQuery(QName);

if nargin > 0,
  Query.Name = QName;
else                        % change the following line to change the query!
  Query.Name = 'BasepairGeometric2';
  Query.Name = 'Sarcin9Mixed';
  Query.Name = 'Sarcin5Geo';
  Query.Name = 'StackedPair'; 
  Query.Name = 'Stack';
  Query.Name = 'Basepair';
  Query.Name = 'Sarcin5Symb';
end

Query.SearchFiles = '1s72';        % default is to search 1s72

switch Query.Name

% ---------------------------------- Searches from FR3D paper by Sarver et al

case 'Sarcin5Geo'
  Query.Description    = 'Sarcin five nucleotide geometric';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};   % all in the 23S
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {'1s72' 'Nonredundant_list'};
  Query.SearchFiles    = {'1s72' 'HighResolution_list'};

case 'Sarcin5Symb'
  Query.Description    = 'Sarcin five nucleotide symbolic';
  Query.Edges{1,2}     = 'tHS';
  Query.Edges{3,5}     = 'cHS';
  Query.Edges{3,4}     = 'tWH';
  Query.MaxDiff(5,3)   = 2;
  Query.MaxDiff(3,1)   = 2;
  Query.MaxDiff(4,2)   = 2;

case 'Sarcin7Mixed'
  Query.Description    = 'Sarcin seven nucleotide mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{3,4}     = 'tWH';
  Query.DiscCutoff     = 0.5;       
  Query.ExcludeOverlap = 0;

case 'Sarcin9Mixed'
  Query.Description    = 'Sarcin nine nucleotide mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692' '2691' '2703' '2690' '2704'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0' '0' '0' '0'};% all in the 23S
  Query.Edges{3,4}     = 'tWH';
  Query.DiscCutoff     = 0.5;
  Query.ExcludeOverlap = 1;
  Query.SearchFiles    = {'1s72' 'Nonredundant_list'};

case 'KinkTurnCentral'
  Query.Description    = 'Kink-turn central base mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '94' '98'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{1,2}     = 'tHS';
  Query.DiscCutoff     = 0.7;  
  Query.ExcludeOverlap = 1;

case 'KinkTurnClosing'
  Query.Description    = 'Kink-turn closing base pair mixed';
  Query.Filename       = '1s72';
  Query.NTList         = {'80' '97' '81' '93' '100' '77'};
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};   % all in the 23S
  Query.Edges{1,2}     = 'tHS';
  Query.DiscCutoff     = 0.9;    
  Query.ExcludeOverlap = 1;

case 'GNRA4NonSeq'
  Query.Description    = 'GRNA hairpin without sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Edges{1,4}     = 'cWW';
  Query.Edges{2,3}     = 'tSH';
  Query.DiscCutoff     = 1;      
  Query.ExcludeOverlap = 1;

case 'GNRA4Seq'
  Query.Description    = 'GRNA hairpin with sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Edges{1,4}     = 'cWW';
  Query.Edges{2,3}     = 'tSH';
  Query.DiscCutoff     = 1;      
  Query.MaxDiff(1,4)   = 6;
  Query.ExcludeOverlap = 1;

case 'GNRA5'
  Query.Description    = 'GRNA hairpin 5 nucleotide';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '805' '807' '808' '809'};
  Query.ChainList      = {'0' '0' '0' '0' '0'}; 
  Query.Edges{1,5}     = 'cWW bif';
  Query.Edges{1,2}     = 's35';
  Query.Edges{3,4}     = 's35';
  Query.DiscCutoff     = 0.8;
  Query.MaxDiff(1,5)   = 6;
  Query.MaxDiff(2,4)   = 4;
  Query.ExcludeOverlap = 1;

% -------------------------------------------- Additional searches

case 'BasepairGeometric'
  Query.Filename   = '1s72';
  Query.NTList     = {'804' '809'};
  Query.DiscCutoff = 0.4;

case 'BasepairGeometric2'
  Query.Filename   = '1u6b';
  Query.NTList     = {'59','85'};
  Query.DiscCutoff = 0.5;
  Query.SearchFiles    = {'Nonredundant_list'};

case 'Basepair'
%  Query.Edges{1,2}     = 'nPB';
  Query.Edges{1,2}      = 'B2P';
%  Query.Diff{1,2}      = '> <1000';
 % Query.Config{1}      = 'syn';
 % Query.SearchFiles    = {'1s72', '1j5e'};

case 'GNRA4NonSeq'
  Query.Description    = 'GRNA hairpin without sequential constraint';
  Query.Filename       = '1s72';
  Query.NTList         = {'804' '809' '805' '808'};
  Query.ChainList      = {'0' '0' '0' '0'}; 
  Query.Mask           = 'NNNN';
  Query.Edges{1,2}     = 'ncWW';
  Query.Edges{3,4}     = 'ntSH';
  Query.DiscCutoff     = 1;         % guaranteed to find all candidates
                                    % with discrepancy below this number
  Query.ExcludeOverlap = 1;

case 'StackedPair'
  Query.Description    = 'Stacked pairs';
  Query.Diff{1,2}      = '> <5';
  Query.Diff{4,3}      = '< <5';
  Query.Edges{1,4}     = 'cWW';
  Query.Edges{2,3}     = 'cWW';
  Query.Edges{1,2}     = 'stack';
  Query.Edges{3,4}     = 'stack';
  Query.Mask           = 'GGNN';

case 'Stack'
  Query.Description    = 'Two stacked bases';
  Query.Mask           = 'AA';
  Query.Edges{1,2}     = 'stack';
  Query.SearchFiles    = {'1s72'};

case 'StackedOncWW'
  Query.Description    = 'What stacks on a cWW?';
  Query.Edges{1,2}     = 'cWW';
  Query.Edges{3,4}     = 'Pair ~tSH ~tHS ~cWW';  % exclude categories we know
  Query.Edges{1,3}     = 'Stack';
  Query.Edges{2,4}     = 'Stack';
  Query.Diff{1,3}      = '=1';
  Query.Diff{2,4}      = '=1';

case 'Sarcin5Quick'
  Query.Description    = 'Sarcin five nucleotide quick search';
  Query.Filename       = '1s72';
  Query.NTList         = {'2694' '2701' '2693' '2702' '2692'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};   % all in the 23S
  Query.DiscCutoff     = 0.3;         % guaranteed to find all candidates
                                      % with discrepancy below this number
case 'cWW-noncWW-cWW'
  Query.Description    = 'Noncanonical pair between canonical';
  Query.Edges{1,6}     = 'cWW';
  Query.Edges{2,5}     = '~cWW';
  Query.Edges{3,4}     = 'cWW';
  Query.MaxDiff(1,2)   = 1;
  Query.MaxDiff(2,3)   = 1;
  Query.MaxDiff(4,5)   = 1;
  Query.MaxDiff(5,6)   = 1;

end

% List of nucleotide mask codes:
% A-A
% C-C
% G-G
% U-U
% A,C-M
% A,G-R
% A,U-W
% C,G-S
% C,U-Y
% G,U-K
% A,C,G-V
% A,C,U-H
% A,G,U-D
% C,G,U-B
% A,C,G,U-N

