% oSpecifyQuery returns the description of a model motif
%
% The variable Query has several fields, most of which are optional:
%   Query.Description    a useful string, can be long
%   Query.Name           a short string, which will become part of a filename
%
% For all searches, these are required:
%   Query.SearchFiles    a cell array of PDB filenames or PDB lists
%   See http://rna.bgsu.edu/rna3dhub/nrlist for current non-redundant lists
%
% For geometric searches, these are required:
%   Query.Filename       a string like 1s72, where the query motif is found
%   Query.NTList         a list of nucleotide numbers
%   Query.ChainList      a list of chain specifications, use if needed
%
% For non-geometric or mixed searches, these are optional:
%   Query.Edges          list of required basepairing or stacking interactions
%                        also specify allowed or disallowed pairs here
%   Query.Mask           a string telling which nucleotides to allow, like 'GNRA'
%     List of nucleotide mask codes:
%     A-A
%     C-C
%     G-G
%     U-U
%     A,C-M
%     A,G-R
%     A,U-W
%     C,G-S
%     C,U-Y
%     G,U-K
%     A,C,G-V
%     A,C,U-H
%     A,G,U-D
%     C,G,U-B
%     A,C,G,U-N
%   Query.AngleWeight    weights to put on the angles (defaults 1)
%   Query.DistanceWeight weights to put on nucleotide distances (defaults 1)
%   Query.DiscCutoff     discrepancy cutoff D_0 (default 0.4)
%   Query.RelCutoff      relaxed cutoff D_1 (default Query.DiscCutoff)
%   Query.MaxDiff        maximum difference between sorted nucleotide indices
%   Query.MinDiff        minimum difference between sorted nucleotide indices
%
%   Query.Geometric      set to 0 to ignore geometry, only use symbolic constraints
%                        Default is 1.
%   Query.ExcludeOverlap set to 1 to eliminate highly redundant motifs with
%     larger discrepancies; often, you get the same candidate with one or two
%     nucleotides different but much higher discrepancy.  Default is 1 when
%     Query.NumNT > 6.

% GG stack
  clear Query
  Query.Name = 'Stacked GG';
  Query.Filename       = '2LX1';
  Query.NTList         = {'6' '8'};
  Query.ChainList      = {'A' 'A'};     %
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {'nrlist_2.134_3.0A_list'};
  Query.SearchFiles    = {'1S72'};
  oFR3DSearch
  return


% AG basepairs

  clear Query              % remove any previous Query parameters
  Query.Name           = 'All AG basepairs';
  Query.Edges{1,2}     = 'AG pair';
  Query.SearchFiles    = {'1S72'};
  oFR3DSearch
  return

% cSS basepairs

  clear Query              % remove any previous Query parameters
  Query.Name           = 'All cSS basepairs';
  Query.Edges{1,2}     = 'csS CC GC UU';
  Query.SearchFiles    = {'1S72','2AW7','3u5h','488d'};
  oFR3DSearch
  return

% AG tHS basepairs
  clear Query              % remove any previous Query parameters
  Query.Name           = 'AG tHS basepairs';
  Query.Edges{1,2}     = 'AG tHS';
  Query.SearchFiles    = {'Ribosome_list'};
  oFR3DSearch
  return

% AA tHH basepairs

  clear Query              % remove any previous Query parameters
  Query.Name           = 'AA tHH basepairs in high resolution structures';
  Query.Edges{1,2}     = 'AA tHH';
  Query.SearchFiles    = {'nrlist_1.52_3.0A_list'};                                       % a locally-stored list
  Query.SearchFiles    = {'http://rna.bgsu.edu/rna3dhub/nrlist/download/1.52/3.0A/csv'};  % download list from this URL
  oFR3DSearch
  return

% Sarcin-Ricin loop core five nucleotides, geometric search

  clear Query              % remove any previous Query parameters
  Query.Name           = 'Sarcin five nucleotide geometric search in ribosomes';
  Query.Filename       = '1s72';
  Query.NTList         = {'2692' '2693' '2694' '2701' '2702'};
  Query.ChainList      = {'0' '0' '0' '0' '0'};     % all in the 23S chain
  Query.DiscCutoff     = 0.5;
  Query.SearchFiles    = {'Ribosome_list'};
  oFR3DSearch
  return

% Sarcin-Ricin loop core five nucleotides, symbolic search

  clear Query          % remove any previous Query parameters
  Query.Name         = 'Sarcin five nucleotide symbolic search in ribosomes';
  Query.Edges{1,2}   = 'cSH';
  Query.Edges{3,4}   = 'tHS';
  Query.Edges{2,5}   = 'tWH';
  Query.Diff{2,1}    = '> =1';
  Query.Diff{3,2}    = '> =1';
  Query.Diff{5,4}    = '> =1';
  Query.SearchFiles  = {'Ribosome_list'};
  oFR3DSearch
  return

% Sarcin-Ricin loop core five nucleotides plus bases beyond the triple, symbolic search

  clear Query              % remove any previous Query parameters
  Query.Name           = 'Sarcin five nucleotide bases beyond triple in ribosomes';
  Query.Edges{5,6}     = 'cSH';
  Query.Edges{7,8}     = 'tHS';
  Query.Edges{6,9}     = 'tWH';
  Query.Diff{6,5}      = '<= 2';
  Query.Diff{7,6}      = '<= 2';
  Query.Diff{9,8}      = '<= 2';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '> =1';
  Query.Diff{4,3}      = '> =1';
  Query.Diff{5,4}      = '> =1';
  Query.Diff{10,9}     = '> =1';
  Query.Diff{11,10}    = '> =1';
  Query.Diff{12,11}    = '> =1';
  Query.Diff{13,12}    = '> =1';
  Query.Diff{14,13}    = '> =1';
  Query.SearchFiles    = {'Ribosome_list'};
  oFR3DSearch
  return

% Four-nucleotide hairpin loop

  clear Query              % remove any previous Query parameters
  Query.Name           = 'Four-nucleotide hairpin loops in ribosomes';
  Query.Edges{1,6}     = 'cWW';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '> =1';
  Query.Diff{4,3}      = '> =1';
  Query.Diff{5,4}      = '> =1';
  Query.Diff{6,5}      = '> =1';
  Query.SearchFiles    = {'Ribosome_list'};
  oFR3DSearch
  return

% Border SS nucleotides making a cWW pair and so forming a hairpin

  clear Query              % remove any previous Query parameters
  Query.Name           = 'cWW Border SS hairpin closing pair in 1S72';
  Query.Edges{1,2}     = 'cWW borderSS';
  Query.Diff{2,1}      = '>';
  Query.SearchFiles    = {'1S72'};
  oFR3DSearch
  return

% Border SS nucleotides making a cWW pair plus more of the hairpin

  clear Query              % remove any previous Query parameters
  Query.Name           = 'cWW Border SS hairpin closing pair plus adjacent nts in 1S72';
  Query.Edges{1,4}     = 'cWW borderSS';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '>';
  Query.Diff{4,3}      = '> =1';
  Query.SearchFiles    = {'1S72'};
  oFR3DSearch
  return

% Internal loops with at least one nucleotide internal to each strand

  clear Query              % remove any previous Query parameters
  Query.Name           = 'Internal loops with at least one nucleotide in 2QBG';
  Query.Edges{1,2}     = 'borderSS';
  Query.Edges{2,3}     = 'cWW';
  Query.Edges{3,4}     = 'borderSS';
  Query.Edges{1,4}     = 'cWW';
  Query.Diff{2,1}      = '>';
  Query.Diff{4,3}      = '>';
  Query.SearchFiles    = {'2QBG'};
  Query.SearchFiles    = {'nrlist_current_4.0A_list'};
  oFR3DSearch
  return

% GNRA hairpin geometric search
  clear Query
  Query.Name           = 'GRNA hairpin geometric search in 1S72';
  Query.Filename       = '1S72';                    % file containing query instance
  Query.NTList         = {'804' '805' '806' '807' '808' '809'}; % nucleotide numbers of query
  Query.ChainList      = {'0' '0' '0' '0' '0' '0'};         % chains of query (optional)
  Query.DiscCutoff     = 0.5;              % limit on geometric discrepancy
  Query.Diff{2,1}      = '>';
  Query.Diff{3,2}      = '>';
  Query.Diff{4,3}      = '>';
  Query.Diff{5,4}      = '>';
  Query.Diff{6,5}      = '>';
  Query.ExcludeOverlap = 1;                % exclude very similar candidates
  Query.SearchFiles    = {'1S72'};
  oFR3DSearch
  return

% End of WC helix, start of some kind of loop
  clear Query              % remove any previous Query parameters
  Query.Name           = 'End of WC helix, start of loop, see nts 3 and 4';
  Query.Edges{1,6}     = 'cWW nested';
  Query.Edges{2,5}     = 'cWW nested';
  Query.Edges{3,4}     = '~cWW';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '> =1';
  Query.Diff{5,4}      = '> =1';
  Query.Diff{6,5}      = '> =1';
  Query.ExcludeOverlap = 1;
  Query.SearchFiles    = {'nrlist_1.52_4.0A_list'};
  oFR3DSearch
  return

% End of WC helix, start of internal loop
  clear Query              % remove any previous Query parameters
  Query.Name           = 'End of WC helix, start of internal loop, nts 2 and 5';
  Query.Edges{1,3}     = 'borderSS';
  Query.Edges{3,4}     = 'cWW';
  Query.Edges{4,6}     = 'borderSS';
  Query.Edges{1,6}     = 'cWW';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '>';
  Query.Diff{5,4}      = '>';
  Query.Diff{6,5}      = '> =1';
  Query.ExcludeOverlap = 1;
  Query.SearchFiles    = {'nrlist_1.46_4.0_list'};
  oFR3DSearch
  return

% End of WC helix, start of junction loop
  clear Query              % remove any previous Query parameters
  Query.Name           = 'End of WC helix, start of junction loop, see nts 2 and 5';
  Query.Edges{1,3}     = 'borderSS';
  Query.Edges{3,4}     = '~cWW';
  Query.Edges{4,6}     = 'borderSS';
  Query.Edges{1,6}     = 'cWW';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '>';
  Query.Diff{5,4}      = '>';
  Query.Diff{6,5}      = '> =1';
  Query.ExcludeOverlap = 1;
  Query.SearchFiles    = {'nrlist_1.46_4.0_list'};
  oFR3DSearch
  return

% End of WC pseudoknot, next bases
  clear Query              % remove any previous Query parameters
  Query.Name           = 'End of WC pseudoknot, start of loop, see nts 3 and 4';
  Query.Edges{2,5}     = 'cWW LR';
  Query.Edges{1,6}     = 'cWW LR';
  Query.Edges{3,4}     = '~cWW';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{3,2}      = '> =1';
  Query.Diff{5,4}      = '> =1';
  Query.Diff{6,5}      = '> =1';
  Query.ExcludeOverlap = 1;
  Query.SearchFiles    = {'nrlist_1.46_4.0_list'};
  oFR3DSearch
  return

% Two helices making a long-range interaction
  clear Query              % remove any previous Query parameters
  Query.Name           = 'Two helices making a long-range interaction';
  Query.Edges{1,2}     = 'stack';
  Query.Edges{3,4}     = 'stack';
  Query.Edges{5,6}     = 'stack';
  Query.Edges{7,8}     = 'stack';
  Query.Edges{1,4}     = 'cWW local';
  Query.Edges{2,3}     = 'cWW local';
  Query.Edges{5,8}     = 'cWW local';
  Query.Edges{6,7}     = 'cWW local';
  Query.Edges{1,8}     = 'LR';
  Query.Diff{2,1}      = '> =1';
  Query.Diff{4,3}      = '> =1';
  Query.Diff{6,5}      = '> =1';
  Query.Diff{8,7}      = '> =1';
  Query.SearchFiles    = {'1S72'};
  Query.SearchFiles    = {'Ribosome_list'};
  oFR3DSearch
  return



% the searches below this line are older and might not work right

switch a
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
  Query.Edges{1,2}      = 'B2P';

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

