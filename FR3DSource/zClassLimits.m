% zClassLimits stores the cutoffs for the computer classification of pairs
% according to displacement, normal vector, and angle of rotation.  As such,
% it is the repository of expert knowledge of pair classifications.

% The main categories have integer values
% The subcategories have decimal parts.
% For decimal parts smaller than 0.5, the same hydrogen bonds will be checked
% as for the main category.
% For decimal parts 0.5 and larger, no hydrogen bonds will be checked

function [ClassLimits, CurrentVersion] = zClassLimits;

CurrentVersion = 8.0;                       % version number of class limits

% 7.3 2011-03-11 Added base-ribose interactions to zInteractionRange
% 7.4 2012-11-21 Added AA 1.1 (it seems it used to be there, was removed?)
% 8.0 2016-12-14 Heavily revised basepair cutoffs, angle in plane, gap restrictions

ClassLimits = zeros(50,12,16);

% For each paircode there is a matrix B which tells the upper and lower
% limits of various variables.
% The first column of B is the category, as in zListCategoryNames
% Some categories are further subdivided with decimal extensions.
% Columns 2 and 3 of B are upper and lower limits on displ(1).
% Columns 4 and 5 of B are upper and lower limits on displ(2).
% Columns 6 and 7 of B are upper and lower limits on displ(3).
% Columns 8 and 9 of B are upper and lower limits on normal(3).
% Columns 10 and 11 of B are limits on Ang (FlipAng or RotAng, depending on
% normal(3)).
% Column 12 is an upper limit on the gap between the edges.  Already somewhat enforced in zAnalyzePair, but this lets you put tighter limits

% AA pairs (paircode 1) ------------------------------------------------------

% =    [  1    5.0  6.5  8.0 10.0 -3.0  3.0 -1.1 -0.7   70  105 5.0];  % cWw
B =    [  1    5.2  7.2  8.3 10.4 -3.0  3.0 -1.1 -0.7   75  105 1.6];  % cWw moved up and right, ang tightened left 2015-06-13
% = [B;[  1.1  6.0  8.2  7.0  8.3 -3.0  3.0 -1.1 -0.7   60   95 5.0]]; % cWw CLZ 2012-11-21
B = [B;[  1.1  6.5  8.2  6.2  8.3 -3.0  3.0 -1.1 -0.7   55   85 1.6]]; % cWw CLZ 2012-11-21; moved right and down, ang moved left 2015-06-15 CLZ
% = [B;[  1.1  6.0  8.0  6.5  8.0 -3.0  3.0 -1.1 -0.7   55   80 5.0]]; % cWW
% = [B;[  2    5.1  8.0  7.7  9.8 -3.0  3.0  0.58 1.1  155  195 5.0]]; % tWW
B = [B;[  2    5.1  8.2  7.7  9.8 -3.0  3.0  0.58 1.1  160  200 1.6]]; % tWW stretched right 2015-06-15 CLZ
% = [B;[  3    5.2  7.2 -1.9 -0.9 -3.5  3.5  0.6  1.1   45   75 5.0]]; % cWH?
% = [B;[  4    3.2  6.6  7.6  9.1 -3.0  3.0 -1.1 -0.7  145  195 5.0]]; % tWH Expanded 2011-07-27 CLZ
B = [B;[  4    3.1  5.8  7.6  9.1 -3.0  3.0 -1.1 -0.7  165  195 1.6]]; % tWH Expanded 2011-07-27 CLZ; shifted left, angle tightened left 2015-06-19 CLZ
% = [B;[  5    3.6  5.5  6.7  9.0 -3.0  3.0 -1.1 -0.6   10   60 5.0]]; % cWS Expanded 2010-11-27 CLZ
B = [B;[  5    3.6  5.7  6.3  9.0 -3.0  3.0 -1.1 -0.6    0   50 1.6]]; % cWS Expanded 2010-11-27 CLZ; expanded right and down, angle expanded left 2015-06-29 CLZ
% = [B;[  5.1  4.6  6.6  5.6  7.3 -3.0  3.0 -1.1 -0.6   10   60 5.0]]; % cWS Added 2010-11-27 CLZ
B = [B;[  5.1  4.6  6.6  5.6  7.1 -3.0  3.0 -1.1 -0.6   -5   25 1.6]]; % cWS Added 2010-11-27 CLZ; top lowered
% = [B;[  5.2  2.0  3.6  7.8  9.0 -3.0  3.0 -1.1 -0.6   40   65 1.6]]; % cWS Added 2011-07-27 CLZ; removed, only one instance 2016-11-11 CLZ
% = [B;[  6    2.2  3.4  8.1  9.1 -3.0  3.0  0.7  1.1  260  -70 5.0]]; % tWS Expanded 2011-07-27 CLZ
B = [B;[  6    1.9  3.4  8.1  9.1 -3.0  3.0  0.7  1.1   75   95 1.6]]; % tWS Expanded 2011-07-27 CLZ; expanded left 2015-06-20 CLZ
B = [B;[  8   -3.2 -0.9  6.4  8.6 -3.0  3.0  0.7  1.1  160  200 1.6]]; % tHh Expanded 2011-07-27 CLZ
% = [B;[  9   -5.5 -2.9  4.6  6.2 -3.0  3.0  0.7  1.1  -25   45 5.0]]; % cHS Expanded 2011-07-27 CLZ
B = [B;[  9   -5.0 -2.9  4.6  6.2 -3.0  3.0  0.7  1.1  -20   10 1.6]]; % cHS Expanded 2011-07-27 CLZ; tightened left, ang shifted lift 2015-06-22 CLZ
% = [B;[ 10   -1.8  0.6  7.1  8.9 -3.0  3.0 -1.1 -0.7   65  110 5.0]]; % tHS Expanded 2011-07-27 CLZ
B = [B;[ 10   -2.0  0.6  7.3  8.9 -3.0  3.0 -1.1 -0.7   65  105 1.6]]; % tHS Expanded 2011-07-27 CLZ
B = [B;[ 11    2.5  5.5 -6.6 -5.0 -3.0  3.1 -1.1 -0.6  265  -55 2.1]]; % cSs 2015-06-23 CLZ
B = [B;[ 12    6.1  8.2  1.1  4.1 -3.0  3.0  0.65 1.1  115  160 2.1]]; % tSs
B = [B;[ 13    6.7  8.6  8.0  9.1 -3.0  3.0 -1.1 -0.7  120  145 1.6]]; % bif

s = size(B);
ClassLimits(1:s(1),1:s(2),1) = B;                      % AA is paircode 1

% AC pairs (paircode 5) ------------------------------------------------------

B =    [  1    4.35 7.0  5.5  8.5 -3.0  3.0 -1.1 -0.7   35   90 1.6];  % cWW
% = [B;[  1.1  6.4  7.8  4.3  6.0 -3.0  3.0 -1.1 -0.7   50   80 5.0]]; % cWW adjusted 2010-07-13 renamed 2010-07-18 removed 2010-11-04 with Neocles
% = [B;[  1.5  8.0  9.9  1.5  4.0 -3.0  3.0 -1.1 -0.7   50   78 5.0]]; % cWW adjusted 2010-07-13 renamed 2010-07-18
B = [B;[  1.5  8.0  9.8  1.5  4.0 -3.0  3.0 -1.1 -0.7   40   75 1.6]]; % cWW tightened right, CLZ 2015-10-11
B = [B;[  2    3.0  4.7  7.8  9.3 -3.0  3.0  0.7  1.1  145  175 1.6]]; % tWW Expanded, angle tightened 2011-07-27 CLZ; Expanded right 2015-08-15 CLZ
% = [B;[ -4   -3.1 -0.2  6.9  9.2 -3.0  3.0 -1.1 -0.6  135  175 5.0]]; % tHW Expanded 2010-07-13, Expanded a bit 2011-07-27
B = [B;[ -4   -3.1 -0.2  6.9  8.7 -3.0  3.0 -1.1 -0.6  140  175 1.6]]; % tHW Top lowered 2015-10-12 CLZ
B = [B;[ -4.1 -4.2 -2.5  5.5  6.8 -3.0  3.0 -1.1 -0.8  150  175 1.6]]; % tHW Shifted down and left, angle tightened right 2015-08-16 CLZ
% = [B;[  5    3.0  6.1  6.5  9.0 -3.0  3.0 -1.1 -0.5   10   54 5.0]]; % cWS
B = [B;[  5    3.0  6.1  6.5  8.7 -3.0  3.0 -1.1 -0.6    0   45 1.6]]; % cWS Top lowered 2015-10-12 CLZ
% = [B;[  5.1  2.0  3.0  8.6  9.4 -3.0  3.0 -1.1 -0.7   50   85 5.0]]; % cWS Adjusted a lot 2011-07-28 CLZ; angle expanded left 2015-11-01 CLZ
B = [B;[  5.1  2.0  3.0  8.6  9.4 -3.0  3.0 -1.1 -0.7   50   65 1.6]]; % cWS Adjusted a lot 2011-07-28 CLZ; angle expanded left 2015-11-01 CLZ
% = [B;[  5.2  5.5  7.2  5.6  7.2 -3.0  3.0 -1.1 -0.6  -50   50 5.0]]; % cWS Added 2011-07-28 CLZ
B = [B;[  5.2  5.5  7.2  5.6  7.2 -3.0  3.0 -1.1 -0.6  -10   20 1.6]]; % cWS Added 2011-07-28 CLZ; angle tightened 2015-11-01 CLZ
% = [B;[ -5    3.5  6.8 -6.1 -4.6 -3.3  3.3 -1.1 -0.6  -60   55 5.0]]; % cSW Expanded 2011-07-28 CLZ
B = [B;[ -5    3.5  6.3 -6.1 -4.6 -3.3  3.3 -1.1 -0.6  -10   20 1.6]]; % cSW Expanded 2011-07-28 CLZ; Tightened right 2015-11-01 CLZ
% = [B;[  6    1.7  3.8  8.0  9.2 -3.0  3.0  0.7  1.1  250  -50 5.0]]; % tWS Adjusted 2011-07-29 CLZ
B = [B;[  6    1.5  2.7  8.2  9.0 -3.0  3.0  0.7  1.1   55   80 1.6]]; % tWS Adjusted 2011-07-29 CLZ; Moved left 2015-11-01 CLZ; Tightened 2016-11-11 CLZ
% = [B;[  6.5 -2.5  1.0  7.8  9.0 -3.0  3.0  0.7  1.1  -60  -40 5.0]]; % tWS Added 2011-07-29 CLZ; These are much closer to cHS, and so belong there
% = [B;[ -6    5.9  7.1 -5.6 -4.4 -3.0  3.0  0.7  1.1   20   70 5.0]]; % tWS Requiring just one h-bond would add 2 instances. 2011-07-29 CLZ
B = [B;[ -6    5.7  7.1 -5.4 -4.0 -3.0  3.0  0.7  1.1  -65  -30 1.6]]; % tSW Moved up and left 2015-11-01 CLZ
% = [B;[  8   -4.7 -2.6  6.4  7.9 -3.0  3.0  0.6  1.1  170  210 5.0]]; % tHH These make no sense. h clashes. 2011-07-29 CLZ
B = [B;[  8   -4.7 -3.0  6.6  8.2 -3.0  3.0  0.6  1.1  150  185 1.6]]; % tHH h clashes. 2011-07-29 CLZ. Moved up 2015-11-01 CLZ
% = [B;[  9   -4.3 -3.1  4.8  6.1 -3.0  3.0  0.6  1.1  -50   50 5.0]]; % cHS Expanded 2011-07-29 CLZ
B = [B;[  9   -4.1 -2.9  4.8  6.1 -3.0  3.0  0.6  1.1  -30   10 1.6]]; % cHS Expanded 2011-07-29 CLZ. Moved right 2016-06-11 CLZ.
% = [B;[  9.5 -3.2 -1.4  6.3  8.2 -3.0  3.0  0.6  1.1  -80  -30 5.0]]; % cHS Added 2011-07-29 CLZ
B = [B;[  9.5 -3.2 -1.4  6.5  8.2 -3.0  3.0  0.6  1.1   10   60 1.6]]; % cHS Added 2011-07-29 CLZ. Tightened bottom 2016-06-11 CLZ.
% = [B;[ -9    2.5  5.5 -6.0 -4.5 -3.0  3.0  0.8  1.1  -40   30 5.0]]; % cSH Expanded 2011-07-29 CLZ
B = [B;[ -9    3.0  5.5 -6.0 -4.5 -3.0  3.0  0.8  1.1  -20   25 1.6]]; % cSH Expanded 2011-07-29 CLZ. Widened angle right 2016-06-11 CLZ.
% = [B;[ 10   -2.5  1.6  7.4  9.5 -3.0  3.0 -1.1 -0.8   60  125 5.0]]; % tHS 2008-05-20 CLZ.  Good 2011-07-29 CLZ
B = [B;[ 10   -2.5  1.0  7.2  9.5 -3.0  3.0 -1.1 -0.8   60  120 1.6]]; % tHS 2008-05-20 CLZ.  Good 2011-07-29 CLZ. Expanded bottom 2016-06-11 CLZ.
B = [B;[-10    6.9  8.4 -4.5 -1.7 -3.0  3.0 -1.1 -0.8   60   85 1.6]]; % tSH Expanded 2011-07-29 CLZ. Expanded top 2016-06-11 CLZ.
B = [B;[ 11.0  2.3  4.7 -7.4 -5.8 -3.0  3.0 -1.1 -0.5  250  -70 3.0]]; % cSs Good 2011-07-29 CLZ. Good 2016-06-11 CLZ.
% = [B;[ 11.1  4.7  7.0 -7.4 -5.3 -3.0  3.0 -1.1 -0.7  235  -65 5.0]]; % cSs with water - but really this is just cSs with a G in the middle! 2010-07-13
B = [B;[-11.0  5.8  7.7 -2.6 -0.3 -3.0  3.0 -1.1 -0.65 260  -55 3.0]]; % csS Expanded 2010-12-01 CLZ.  OK, but some are in different planes 2011-07-29 CLZ
% = [B;[-11.1  3.9  5.8 -3.6 -2.1 -3.0  3.0 -1.1 -0.65 225  -80 5.0]]; % csS but really out of plane and not a basepair
% = [B;[-11.1  5.8  7.3 -3.3 -2.6 -3.0  3.0 -1.1 -0.65 250  -50 5.0]]; % csS Jesse 10_11_07, modified by CLZ 10-16-07. Only 4 instances 2011-07-29 CLZ
B = [B;[-11.1  6.4  7.3 -3.3 -2.6 -3.0  3.0 -1.1 -0.75 260  -80 1.6]]; % csS Jesse 10_11_07, modified by CLZ 10-16-07. Only 4 instances 2011-07-29 CLZ. Tightened left 2016-06-11 CLZ.
B = [B;[ 12    6.5  8.7  1.0  5.0 -3.0  3.0  0.55 1.1  100  160 2.0]]; % tSs Expanded slightly, jumbled, this is what near category is for 2011-07-29 CLZ. Good 2016-06-11 CLZ.
B = [B;[ 13    7.0 10.5  4.0  7.8 -3.0  3.0 -1.1 -0.7  100  150 1.6]]; % bif Greatly expanded 2011-07-29 CLZ. Angle narrowed left 2016-06-11 CLZ.
B = [B;[-13    1.5  2.8  8.3  9.8 -3.0  3.0 -1.1 -0.7   85  135 1.6]]; % bif Expanded 2011-07-29 CLZ. Angle expanded left 2016-06-11 CLZ.
B = [B;[-13.1  3.0  4.7  7.9  9.3 -3.0  3.0 -1.1 -0.7   45  100 1.6]]; % bif Expanded 2011-07-29 CLZ. Good 2016-06-11 CLZ.

s = size(B);
ClassLimits(1:s(1),1:s(2),5) = B;                      % AC is paircode 5

% CC pairs (paircode 6) ------------------------------------------------------

% =    [  1    6.3  7.2 -0.7  1.3 -2.9  2.9 -1.1 -0.7   30   65 5.0];  % cWw      Nearly all ncWW are cross-strand in helix. Meaningless!
B =    [  1    6.6  7.2 -0.8  1.6 -2.9  2.9 -1.1 -0.7   30   60 1.5];  % cWw Shifted right, expanded vertically, angle tightened 2016-06-11
% = [B;[  2    4.3  5.8  4.3  6.5 -3.0  3.0  0.7  1.1  160  190 5.0]]; % tWW - JS_7_4_08
B = [B;[  2    4.3  5.8  4.2  5.7 -3.0  3.0  0.7  1.1  170  200 1.6]]; % tWW - JS_7_4_08  Shifted down 2016-06-12.
B = [B;[  2.1  3.2  5.1  6.3  7.8 -3.0  3.0  0.7  1.1  165  190 1.6]]; % tWW New subcategory 2016-06-13 CLZ
B = [B;[  3    6.8  8.2  1.0  3.6 -3.0  3.0  0.7  1.1  -90  -70 1.6]]; % cWH
B = [B;[  3.5  4.9  6.4  0.2  1.9 -3.0  3.0  0.7  1.1  -40  -10 1.6]]; % cWH New subcategory 2016-06-19 CLZ
% = [B;[  4    5.6  7.6  2.0  5.0 -3.0  3.0 -1.1 -0.7  115  165 5.0]]; % tWH
B = [B;[  4    5.8  7.4  3.0  6.1 -3.0  3.0 -1.1 -0.7  140  165 1.6]]; % tWH Tightened bottom, left, right 2016-06-19 CLZ
B = [B;[  5    3.6  5.2  5.1  6.8 -3.0  3.0 -1.1 -0.7  -15   15 1.6]]; % cWS Tightened right, angle widened 2016-06-19 CLZ
% = [B;[  5.1  2.2  4.8  7.0  8.6 -3.0  3.0 -1.1 -0.7   70   80 5.0]]; % cWS? Merged with bif category 2010-07-13
% = [B;[  6    2.4  3.5  7.1  8.4 -3.0  3.0  0.5  1.1  240  -80 5.0]]; % tWS
B = [B;[  6    2.4  3.5  7.1  8.4 -3.0  3.0  0.6  1.1   85  110 1.6]]; % tWS Tightened normal vector 2016-06-19 CLZ
B = [B;[  6.5  3.5  5.0  6.8  8.4 -3.0  3.0  0.6  1.1   95  130 1.6]]; % tWS  New subcategory 2016-06-19 CLZ
% = [B;[  8   -4.7 -2.6  6.4  7.5 -3.0  3.0  0.6  1.1  170  210 5.0]]; % tHh Removed 12-14-2007 JS
B = [B;[  9   -6.0 -4.0  3.0  5.4 -3.0  3.0  0.6  1.1  -25   25 1.6]]; % cHS OK with just one hydrogen, 2010-12-21 CLZ
B = [B;[  9.6 -7.2 -6.2  1.5  3.0 -3.0  3.0  0.6  1.1  -30    0 1.6]]; % cHS New subcategory 2016-06-19 CLZ
B = [B;[  9.5 -4.5 -2.5  6.4  8.5 -3.0  3.0  0.6  1.1   10   75 1.6]]; % cHS CLZ 2010-04-29 Tightened 2010-12-21
B = [B;[ 10   -2.5 -0.4  6.5  8.5 -3.0  3.0 -1.1 -0.7   60  100 1.6]]; % tHS CLZ 2010-10-27 but hydrogen bonds will disqualify some; Shifted right 2016-06-19 CLZ
B = [B;[ 11    2.6  4.4 -7.5 -5.2 -3.0  3.0 -1.1 -0.7  260  -65 1.6]]; % cSs CLZ 2008-04-22. Moved up and right, angle moved right 2016-06-19 CLZ
% = [B;[ 13    7.2  9.0 -0.8  3.9 -3.0  3.0 -1.1 -0.7   55   95 5.0]]; % bif updated 2010-07-13 to include cWSa 5.1.
B = [B;[ 13    7.2  9.0 -0.8  3.9 -3.0  3.0 -1.1 -0.7   40   95 1.6]]; % bif updated 2010-07-13 to include cWSa 5.1. Angle tightened left 2016-06-19 CLZ
% = [B;[ 13.1  6.9  7.9  4.3  5.4 -3.0  3.0 -1.1 -0.7   60   80 5.0]]; % based on a strange feature of 1QCU

s = size(B);
ClassLimits(1:s(1),1:s(2),6) = B;                      % CC is paircode 6

% GC pairs (paircode 7) ------------------------------------------------------
% Note:  Before 2016-11-11, the signs of the numbers for the categories were reversed, -5 instead of 5, and special code made it work out right.
% After 2016-11-11, this should not be necessary, but it may take some time to work out all of the bugs.

% =    [  1    5.6  9.0  2.8  6.9 -3.0  3.0 -1.1 -0.7   45   90 5.0];  % cWW
B =    [  1    5.6  9.0  3.1  6.9 -3.0  3.0 -1.1 -0.7   45   90 1.6];  % cWW Fine tuning 2016-06-19 CLZ
% = [B;[  1.1  5.4  6.3  6.6  7.6 -3.0  3.0 -1.1 -0.7   70  110 5.0]]; % cWW
B = [B;[  1.1  4.5  6.3  6.5  7.7 -3.0  3.0 -1.1 -0.7   70  100 1.6]]; % cWW Expanded up and right, angle tightened 2016-07-15 CLZ
% = [B;[  2    6.0  8.3  4.2  6.4 -3.0  3.0  0.6  1.1  190  220 5.0]]; % tWW CLZ 2010-07-15
B = [B;[  2    6.0  8.3  4.2  6.7 -3.0  3.0  0.6  1.1  145  170 1.6]]; % tWW CLZ 2010-07-15
% = [B;[ -3   -4.8 -2.1  4.5  8.0 -3.0  3.0  0.7  1.1  245  -50 5.0]]; % cHW
B = [B;[ -3   -5.1 -3.5  4.0  6.4 -3.0  3.0  0.7  1.1   55   90 1.6]]; % cHW This is the LSW category 2016-07-16 CLZ (concern about 3 versus -3)
B = [B;[ -3.1 -4.0 -2.0  6.4  8.0 -3.0  3.0  0.7  1.1   80  110 1.6]]; % cHW This is a bit more open 2016-07-16 CLZ
% = [B;[ -4   -1.9 -0.7  7.6  8.8 -3.0  3.0 -1.1 -0.7  135  175 5.0]]; % tHW
B = [B;[ -4   -1.9 -0.7  7.8  8.7 -3.0  3.0 -1.1 -0.7  145  175 1.6]]; % tHW Moved up 2016-07-16 CLZ
% = [B;[ -5    4.2  5.2 -6.3 -5.3 -3.0  3.0 -1.1 -0.6   30   60 5.0]]; % cSW
B = [B;[ -5    4.1  5.4 -6.4 -5.3 -3.0  3.0 -1.1 -0.6  -10   10 1.6]]; % cSW Expanded 2016-07-16 CLZ
% = [B;[  5    6.5  7.6  4.0  5.4 -3.0  3.0 -1.1 -0.7   15   40 5.0]]; % cWS
B = [B;[  5    6.6  8.3  3.8  5.2 -3.0  3.0 -1.1 -0.7    5   35 1.6]]; % cWS Moved down and right 2016-07-25 CLZ
B = [B;[ -6    6.2  8.1 -4.0 -2.0 -3.0  3.0  0.7  1.1  255  -65 1.6]]; % tSW Expanded 2010-12-01 CLZ
% = [B;[  6    4.7  6.2  5.9  7.3 -3.0  3.0  0.7  1.1  235  -80 5.0]]; % tWS
B = [B;[  6    4.7  6.4  5.9  7.3 -3.0  3.0  0.7  1.1   85  125 1.6]]; % tWS Expanded right 2016-07-25 CLZ
B = [B;[  7   -3.8 -2.6  4.7  5.9 -3.0  3.0 -1.1 -0.7  -30  -10 1.6]]; % cHH
B = [B;[  7.5 -5.0 -3.6  4.1  6.8 -3.0  3.0 -1.1 -0.7  -90  -35 1.6]]; % cHH New category to include nearby basepairs 2016-11-11 CLZ
B = [B;[  8   -4.8 -3.0  6.5  7.7 -3.0  3.0  0.7  1.1  150  185 1.6]]; % tHH Tightened CLZ 2010-11-05
% = [B;[  8.1 -2.2  1.1  8.4 11.2 -3.0  3.3  0.7  1.1  160  210 5.0]]; % tHH Expanded CLZ 2010-07-15 A little more 2010-12-06
B = [B;[  8.1 -2.2  1.1  8.2 10.3 -3.0  3.3  0.7  1.1  180  200 1.6]]; % tHH Moved down, angle tightened left 2016-07-25 CLZ
B = [B;[  9    3.7  5.2 -7.7 -6.6 -3.0  3.0  0.7  1.1  -25   -5 1.6]]; % cHS CLZ 2007-12-14
B = [B;[  9.1  2.4  4.0 -6.6 -5.5 -3.0  3.0  0.7  1.1   -5   25 1.6]]; % cHS Instances found, expanded 2016-11-11 CLZ
B = [B;[ 11    5.5  7.4 -3.4 -1.9 -3.0  3.0 -1.1 -0.7  250  -60 2.5]]; % cSs Shifted 2010-11-24 CLZ; hard to expand to add instances 2016-07-26
% = [B;[ 11.1  5.7  6.1 -4.0 -3.2 -3.0  3.0 -1.1 -0.7  -50  -20 5.0]]; % cSs Removed because no instances 2016-08-02 CLZ
% = [B;[ 11.2  4.3  5.3 -4.3 -3.3 -3.0  3.0 -1.1 -0.7  -90  -50 5.0]]; % cSs Added 2010-11-24 CLZ for an example from Neocles
B = [B;[ 11.1  4.5  5.5 -4.0 -3.0 -3.0  3.0 -1.1 -0.7  -90  -50 1.6]]; % cSs Moved up and right 2016-08-02 CLZ; no instances 2016-11-11 CLZ
B = [B;[ 11.5  4.0  5.4 -3.7 -2.7 -3.0  3.0 -1.1 -0.7  235  270 1.6]]; % cSs Added 2010-11-26 CLZ
% = [B;[-11    2.6  3.9 -7.1 -5.6 -3.0  3.0 -1.1 -0.7  250  -70 5.0]]; % cSS
B = [B;[-11    1.7  3.4 -8.7 -6.3 -3.0  3.0 -1.1 -0.7  230  270 2.6]]; % csS CLZ 2007-12-13
% = [B;[-12    7.8  9.0 -1.9 -0.5 -3.0  3.0  0.7  1.1  140  200 5.0]]; % tsS
B = [B;[-12    7.7  9.0 -2.1 -0.5 -3.0  3.0  0.7  1.1  165  210 1.6]]; % tsS Expanded down and left 2016-08-05 CLZ

s = size(B);
ClassLimits(1:s(1),1:s(2),7) = B;

% AG pairs (paircode 9) ------------------------------------------------------

% =    [  1    7.0  9.0  5.6  8.4 -3.5  3.5 -1.1 -0.65  65  100 5.0];  % cWW
B =    [  1    7.0  9.0  5.6  8.6 -3.5  3.5 -1.1 -0.65  65  100 1.5];  % cWW Gap reduced 2016-08-05 CLZ
% = [B;[  1.1  8.6  9.6  3.2  4.6 -3.0  3.0 -1.1 -0.7   50   75 5.0]]; % cWW
B = [B;[  1.1  8.8  9.8  3.1  5.1 -3.0  3.0 -1.1 -0.7   50   70 1.6]]; % cWW Adjusted 2016-10-01 CLZ
% = [B;[  1.2  5.8  7.0  7.6  8.6 -3.0  3.0 -1.1 -0.5   90  125 5.0]]; % cWW
B = [B;[  1.2  5.5  7.0  7.8  9.6 -3.0  3.0 -1.1 -0.5   80  115 1.6]]; % cWW Adjusted 2016-10-01 CLZ
% = [B;[  3    6.7  7.8  4.4  5.4 -3.0  3.0  0.7  1.1   70  100 5.0]]; % cWH OK 2010-07-21
B = [B;[  3    6.7  7.8  3.8  5.4 -3.0  3.0  0.7  1.1  260  -75 1.6]]; % cWH OK 2010-07-21, Expanded down, angle, gap tightened 2016-10-01 CLZ
% = [B;[  3.1  4.4  5.4  7.0  8.3 -3.0  3.04 0.7  1.1   75  110 5.0]]; % cWH
B = [B;[  3.5  4.4  5.4  7.6  8.8 -3.0  3.04 0.7  1.1  225  265 1.6]]; % cWH Moved up, angle moved up, renamed cWHe 2016-10-01 CLZ
% = [B;[ -3   -6.0 -4.9  6.5  7.8 -3.0  3.0  0.7  1.1  250  -80 5.0]]; % cHW
B = [B;[ -3   -6.0 -4.9  6.3  7.8 -3.0  3.0  0.7  1.1   80  110 1.6]]; % cHW Expanded down 2016-10-01 CLZ
% = [B;[  4    4.3  5.4  7.2  8.5 -3.0 3.02 -1.1 -0.7  155  185 5.0]]; % tWH
B = [B;[  4    3.6  5.2  7.6  8.8 -3.0 3.02 -1.1 -0.7  165  190 1.6]]; % tWH Moved up and left 2016-10-01 CLZ
% = [B;[  5    3.7  4.7  6.9  8.7 -3.0  3.0 -1.1 -0.5   15   65 5.0]]; % cWS
B = [B;[  5    3.0  5.4  6.9  8.7 -3.0  3.0 -1.1 -0.5    5   45 1.6]]; % cWS Expanded left and right, angle tightened 2016-10-01 CLZ
% = [B;[ -5    6.6  8.6 -4.8 -2.8 -3.0  3.0 -1.1 -0.5  -90   -5 5.0]]; % cSW Jesse 10_11_07 Expanded angle right 2016-10-21 CLZ
B = [B;[ -5    7.0  9.0 -7.0 -3.5 -3.0  3.0 -1.1 -0.5  -55   15 1.6]]; % cSW Jesse 10_11_07 Expanded angle right 2016-10-21 CLZ
B = [B;[  6    0.0  4.0  7.7  9.3 -3.0  3.0  0.6  1.1   50   90 5.0]]; % tWS Expanded 2010-11-27 CLZ
% = [B;[  6.1  6.3  7.8  6.4  8.5 -3.0  3.0  0.7  1.1  220  265 5.0]]; % tWS CLZ removed 2010-07-15; doesn't contribute enough, confuses triples
% = [B;[  7   -4.6 -3.6  7.8  8.8 -3.0  3.0 -1.1 -0.7  180  225 5.0]]; % cHH
B = [B;[  7   -4.8 -2.7  8.2 10.0 -3.0  3.0 -1.1 -0.7  185  220 1.4]]; % cHH Greatly expanded 2016-10-21 CLZ
% = [B;[  8   -3.3 -1.5  8.5 10.2 -3.0  3.02 0.7  1.1  195  245 5.0]]; % tHH
B = [B;[  8   -2.6 -0.3  8.5 10.7 -3.0  3.02 0.7  1.1  115  175 1.6]]; % tHH Moved up and right 2016-10-21 CLZ
% = [B;[  9   -6.6 -4.9  1.5  3.6 -3.0  3.0  0.7  1.1   10   70 5.0]]; % cHS
B = [B;[  9   -7.0 -4.9  0.3  4.0 -3.0  3.0  0.7  1.1  -40   -5 1.6]]; % cHS Expanded 2016-10-21 CLZ
% = [B;[  9.1 -5.2 -2.2  4.5 6.3  -3.0  3.0  0.7  1.1    0   70 5.0]]; % cHS CLZ 2010-10-27; no instances 2016-11-11 CLZ
% = [B;[  9.2 -2.6 -0.5  7.0 8.8  -3.0  3.0  0.7  1.1  -60   40 5.0]]; % cHS Jesse 4_21_08; no instances 2016-11-11 CLZ
% = [B;[ -9    3.5  7.0 -5.1 -3.1 -3.0  3.0  0.7  1.1    0   60 5.0]]; % cSH Jesse 4_18_08. Only one really good case 2011-08-08 CLZ
B = [B;[ -9    4.0  8.0 -5.5 -3.5 -3.0  3.0  0.7  1.1  -45    0 1.6]]; % cSH Jesse 4_18_08. Only one really good case 2011-08-08 CLZ. Moved right 2016-10-21 CLZ
B = [B;[ -9.1  7.7 11.2 -3.8 -0.2 -3.0  3.0  0.7  1.1  -80  -35 1.6]]; % cSH CLZ 2010-10-27
% = [B;[ 10   -2.0  0.7  7.1  9.0 -3.0  3.0 -1.1 -0.6   55  100 5.0]]; % tHS
B = [B;[ 10   -2.0  0.7  7.1  9.0 -3.0  3.0 -1.1 -0.6   50   95 1.6]]; % tHS Reduced gap 2016-10-26 CLZ
% = [B;[ 10.1 -3.3 -2.0  7.1  9.0 -3.0  3.0 -1.1 -0.6   55  100 5.0]]; % tHS
B = [B;[ 10.1 -3.5 -2.0  7.1  9.0 -3.0  3.0 -1.1 -0.6   75  110 1.6]]; % tHS Expanded left, angle moved right 2016-10-30 CLZ
B = [B;[ 10.2  0.3  2.2  8.0  9.4 -3.0  3.0 -1.1 -0.7   40   75 1.6]]; % tHS Shifted down, angle tightened 2016-11-01 CLZ
B = [B;[ 11    2.0  7.0 -6.8 -4.8 -3.0  3.0 -1.1 -0.5  240  -55 2.5]]; % cSs Greatly expanded 2010-11-25 CLZ
B = [B;[-11    5.2  7.5 -2.5  0.5 -3.0  3.0 -1.1 -0.6  230  -40 2.7]]; % csS
B = [B;[ 12    5.4  8.3 -0.6  3.5 -4.2  3.5  0.5  1.1  115  170 1.7]]; % tSs
B = [B;[ 12.1  7.3  8.5  2.5  3.8 -4.2  3.5  0.5  1.1  120  140 1.6]]; % tSs New 2016-11-02 CLZ
% = [B;[-12    7.4  8.7 -1.3 -0.1 -3.0  3.0  0.7  1.1  170  195 5.0]]; % tSS
B = [B;[-12    7.4  9.0 -2.0  0.4 -3.0  3.0  0.7  1.1  165  190 1.7]]; % tSS  Expanded up, down, right 2016-11-02 CLZ
B = [B;[ 14    3.5  6.5 10.0 11.5 -3.0  3.0 -1.1 -0.7  110  140 1.6]]; % water-inserted AG CLZ 2010-09-13; removed huge gap 2016-11-02 CLZ
% = [B;[ 15    5.8  7.7 -5.2 -2.3 -2.5  2.0  0.3  0.9  170  200 5.0]]; % Rib CLZ
B = [B;[ 15    5.8  7.7 -5.2 -2.3 -2.5  2.0  0.3  0.9  165  205 2.4]]; % Rib CLZ  Tightened angle right 2016-11-02 CLZ

s = size(B);
ClassLimits(1:s(1),1:s(2),9) = B;                      % AG is paircode 9

% GG pairs (paircode 11) -----------------------------------------------------

% =    [  1    5.4  6.5  7.7  9.6 -3.6  3.6 -1.1 -0.5   70  105 5.0];  % cWW These only exist when modeled wrong!
% =    [  2    5.2  7.2  7.8  9.8 -3.0  3.0  0.7  1.1  160  210 5.0];  % tWW Jesse 10_11_07
B =    [  2    5.2  7.2  7.7  9.8 -3.0  3.0  0.7  1.1  165  190 1.6];  % tWW Jesse 10_11_07 Tightened right, angle tightened right 2016-11-02 CLZ
B = [B;[  3    7.9 10.0  1.4  4.8 -3.0  3.0  0.7  1.1  255  -70 1.6]]; % cWH
% = [B;[  4    4.6  7.2  5.4  7.4 -3.0  3.0 -1.1 -0.6  155  200 5.0]]; % tHW Expanded 2010-11-26 CLZ
B = [B;[  4    4.4  7.6  5.7  7.4 -3.0  3.0 -1.1 -0.6  160  195 1.6]]; % tHW Expanded 2010-11-26 CLZ  Expanded up and left 2016-11-02 CLZ
% = [B;[  5    5.4  8.1  5.9  7.5 -3.0  3.0 -1.1 -0.7    0   45 5.0]]; % cWS CLZ 2008-07-03
B = [B;[  5    5.2  8.8  5.2  7.5 -3.0  3.0 -1.1 -0.7  -10   20 1.6]]; % cWS CLZ 2008-07-03  Expanded 2016-11-02 CLZ
% = [B;[  7   -2.1 -0.9  5.5  6.7 -3.0  3.0 -1.1 -0.7  -45    0 5.0]]; % cHh
B = [B;[  7   -3.1 -1.3  5.4  6.3 -3.0  3.0 -1.1 -0.7  -20    0 1.8]]; % cHh Moved down and left, angle tightened 2016-11-05 CLZ
B = [B;[  7.1 -2.5 -1.0  7.0  8.0 -3.0  3.0 -1.1 -0.7  -65  -25 1.6]]; % cHh New group, not sequentially adjacent 2016-11-05 CLZ
B = [B;[  8   -0.8  0.4  7.5  8.7 -3.0  3.0  0.7  1.1  235  265 1.6]]; % tHh
% = [B;[  9   -6.3 -4.2  4.8  7.0 -3.0  3.0  0.7  1.1  -45    0 5.0]]; % cHS
B = [B;[  9   -6.3 -4.2  3.5  7.5 -3.0  3.0  0.7  1.1  -25   20 1.6]]; % cHS Expanded up, down, angle moved right 2016-11-05 CLZ
B = [B;[  9.5 -5.5 -3.0  6.8 10.0 -3.0  3.0  0.6  1.1    0   60 1.6]]; % cHS CLZ 2010-11-23 has different hydrogen bond Tightened 2010-12-20
% = [B;[ 10   -1.6  1.3  8.8 10.4 -3.0  3.0 -1.1 -0.7   70  105 5.0]]; % tHS
B = [B;[ 10   -1.6  1.3  8.8 10.8 -3.0  3.0 -1.1 -0.7   65  105 1.6]]; % tHS Expanded up 2016-11-05 CLZ
% = [B;[ 10.1  1.3  2.0  9.1  9.9 -3.0  3.0 -1.1 -0.7   55   75 5.0]]; % tHS
B = [B;[ 10.1  1.3  4.5  9.1 10.6 -3.0  3.0 -1.1 -0.7   35   80 1.2]]; % tHS Expanded right, up, angle expanded, tight gap imposed 2016-11-05 CLZ
% = [B;[ 11    1.0  3.4 -6.9 -4.9 -3.0  3.0 -1.1 -0.7  230  265 5.0]]; % cSs
B = [B;[ 11    0.0  3.0 -8.0 -5.0 -3.0  3.0 -1.1 -0.7  235  265 2.6]]; % cSs Moved down and left 2016-11-05 CLZ
B = [B;[ 12    6.5  8.5 -3.0 -0.2 -3.0  3.0  0.6  1.1  165  205 5.0]]; % tSs Reduced h-bond requirement 2010-12-01 CLZ
% = [B;[ 13    6.6  8.0  7.3  8.6 -3.0  3.0 -1.1 -0.7  125  150 5.0]]; % bif
B = [B;[ 13    6.6  7.8  7.3  8.4 -3.0  3.0 -1.1 -0.7  125  150 1.6]]; % bif  Tightened top, right 2016-11-05 CLZ.  Many near bif, cross strand helix poorly modeled
B = [B;[ 15    5.8  8.0 -6.0 -3.8 -3.5  0.5  0.2  0.8  180  210 2.5]]; % rib CLZ

s = size(B);
ClassLimits(1:s(1),1:s(2),11) = B;                      % GG is paircode 11

% AU pairs (paircode 13) ------------------------------------------------------

B =    [  1    5.6  8.5  3.2  6.8 -3.0  3.0 -1.1 -0.8   45   90 1.6];  % cWW
B = [B;[  2    4.8  6.4  6.0  7.8 -3.0  3.0  0.7  1.1  145  180 1.6]]; % tWW
% = [B;[  3    6.0  8.0  6.0  8.0 -3.0  3.0  0.7  1.1   80  120 5.0]]; % cWH based on a modeled pair
B = [B;[  3    6.5  8.5  4.5  7.2 -3.0  3.0  0.7  1.1  265  -70 1.6]]; % cWH Moved down, angle moved left, actual instances found 2016-11-08 CLZ
B = [B;[ -3   -6.0 -4.0  3.0  5.5 -3.0  3.0  0.7  1.1   50   80 1.6]]; % cHW
% = [B;[ -3.1 -6.5 -5.2  6.5  7.5 -3.0  3.0  0.7  1.1  -55  -40 1.6]]; % cHW No instances found 2016-11-09 CLZ
B = [B;[ -4   -5.0 -1.9  5.0  7.8 -3.0  3.0 -1.1 -0.7  145  185 1.6]]; % tHW
B = [B;[  5    3.3  6.8  6.0  9.0 -3.0  3.0 -1.1 -0.6    0   50 2.0]]; % cWS Expanded 2010-11-25; needs bigger gap than usual
% = [B;[ -5    3.8  5.2 -4.8 -3.3 -3.0  3.0 -1.1 -0.7  -55  -15 5.0]]; % cSW
B = [B;[ -5    3.8  5.3 -4.9 -3.2 -3.0  3.0 -1.1 -0.7  -45   -5 1.6]]; % cSW Expanded 2016-11-09 CLZ
B = [B;[  6   -0.7  2.7  7.8  9.6 -3.0  3.0  0.6  1.1   40   80 1.6]]; % tWS Exchanged with 6.1 CLZ 12-14-2007 - Jesse_4_18_08
% = [B;[ -6    5.8  7.0 -3.2 -1.8 -3.0  3.0  0.7  1.1   90  120 5.0]]; % tSW
B = [B;[ -6    5.8  7.1 -3.6 -1.8 -3.0  3.0  0.7  1.1  250  -80 1.6]]; % tSW Expanded 2016-11-09 CLZ
% = [B;[  8   -3.2 -1.0  8.3 10.5 -3.0  3.0  0.8  1.1  180  210 5.0]]; % tHh
B = [B;[  8   -3.3 -0.6  8.2 10.5 -3.0  3.0  0.7  1.1  145  180 1.6]]; % tHh Expanded, allowed lower normal 2016-11-09 CLZ
% = [B;[  8.1 -4.0 -3.0  5.5  6.5 -3.0  3.0  0.7  1.1  163  183 5.0]]; % tHH No instances 2016-11-09 CLZ
% = [B;[  8.2 -0.8  1.8  9.5 11.0 -3.0  3.0  0.8  1.1  180  220 5.0]]; % tHH No instances 2016-11-09 CLZ
% = [B;[  9.1 -4.2 -3.5  4.5  5.8 -3.0  3.0  0.6  1.1   10   35 5.0]]; % cHS
B = [B;[  9   -4.5 -3.0  4.3  6.0 -3.0  3.0  0.6  1.1  -15   15 1.6]]; % cHS Greatly expanded, renamed as 9.0 2016-11-10-CLZ
% = [B;[  9   -3.2 -1.8  6.2  7.6 -3.0  3.0  0.6  1.1  -60  -35 5.0]]; % cHS
B = [B;[  9.1 -3.5 -1.8  6.2  7.6 -3.0  3.0  0.6  1.1   10   50 1.6]]; % cHS Expanded left, angle expanded right, renamed as 9.1 2016-11-09 CLZ
B = [B;[ -9    3.2  5.6 -5.4 -4.2 -3.0  3.0  0.7  1.1   -5   25 1.6]]; % cHS New, to capture recurrent adjacent basepairing 2016-11-10 CLZ
% = [B;[ -9.1  6.6  8.0 -4.5 -3.5 -3.0  3.0  0.7  1.1   40   80 5.0]]; % cHS Move to category -9.1 2016-11-10 CLZ; no instances
% = [B;[ -4.1 -2.1  0.0  7.7  9.5 -3.0  3.0 -1.1 -0.8   98  155 5.0]]; % tHW Removed 2010-11-06 CLZ too much overlap with 10
% = [B;[ 10    0.5  3.5  8.4 10.3 -3.0  3.0 -1.1 -0.8   75  110 1.6]]; % tHS 2010-11-06 CLZ
% = [B;[ 10.1 -2.5  1.0  7.2  9.5 -3.0  3.0 -1.1 -0.8   60  130 1.6]]; % tHS New 2016-11-10 CLZ. Usually participates in base triple with U
B = [B;[ 10   -2.5  1.0  7.2  9.5 -3.0  3.0 -1.1 -0.8   60  110 1.6]]; % tHS 2010-11-06 CLZ Angle tightened right 2016-11-27 CLZ
B = [B;[ 10.5 -6.5 -4.4  1.0  3.3 -3.0  3.0 -1.1 -0.7  115  165 1.6]]; % tHS renamed 2010-07-18, Top lowered 2016-11-10 CLZ. Far Hoogsteen side of the A.
B = [B;[ 10.6  0.5  3.5  8.4 10.3 -3.0  3.0 -1.1 -0.8   75  105 1.6]]; % tHS New 2016-11-10 CLZ. Usually participates in base triple with U. Different h-bonds.
B = [B;[-10    9.3 11.0 -0.7  1.6 -3.0  3.0 -1.1 -0.7   55  110 1.6]]; % tHS
% = [B;[ 11    3.8  5.3 -6.2 -3.4 -3.0  3.0 -1.1 -0.7  260  -55 5.0]]; % cSs
B = [B;[ 11    3.6  5.6 -6.4 -4.1 -3.0  3.0 -1.1 -0.7  260  -40 3.6]]; % cSs Moved down, widened, angle tightened left 2016-11-10 CLZ
B = [B;[-11    5.7  7.3 -2.5  0.0 -3.5  3.9 -1.1 -0.6  260  -50 3.6]]; % csS Some large base gaps but still sugar-sugar contacts
B = [B;[-11.1  5.3  7.2 -3.2 -2.5 -3.5  3.9 -1.1 -0.6  250  270 1.6]]; % csS New, at the bottom of 11 2016-11-10 CLZ
B = [B;[ 12    6.5  9.2  0.5  4.5 -3.0  3.0  0.7  1.1  105  160 5.0]]; % tSs Expanded 2010-11-24 CLZ

s = size(B);
ClassLimits(1:s(1),1:s(2),13) = B;

% CU pairs (paircode 14) -----------------------------------------------------

B =    [  1    5.6  7.1  5.5  7.2 -3.0  3.0 -1.1 -0.7   80  105 1.6];  % cWW Angle widened right 2016-11-11 CLZ
B = [B;[  1.1  5.0  6.5  6.8  8.4 -3.0  3.0 -1.1 -0.7   100 135 1.4]]; % cWW Angle widened right 2016-11-11 CLZ
B = [B;[  1.2  5.7  6.8  2.3  4.0 -3.0  3.0 -1.1 -0.7   35   85 1.3]]; % cWW Gap tightened a lot; too much variability in these 2016-11-11 CLZ
B = [B;[  1.3  6.0  7.6  3.6  5.7 -3.0  3.0 -1.1 -0.7   70  100 1.4]]; % cWW New subcategory 2016-11-11 CLZ
B = [B;[  1.4  4.0  5.5  7.6  9.0 -3.0  3.0 -1.1 -0.7  140  160 1.4]]; % cWW New subcategory 2016-11-11 CLZ
B = [B;[  2    4.5  5.7  4.3  5.3 -3.0  3.0  0.7  1.1  160  190 1.6]]; % tWW Moved down, left 2016-11-11 CLZ
% = [B;[  3    5.4  8.1  2.0  4.0 -3.0  3.0  0.7  1.1   45   95 1.6]]; % cWH CLZ 2007-12-13
B = [B;[  3    5.8  7.7  2.2  4.8 -3.0  3.0  0.7  1.1  250  -50 1.6]]; % cWH CLZ 2007-12-13; expanded up, angle expanded right 2016-11-11 CLZ
B = [B;[  5    3.7  5.5  5.0  6.6 -3.0  3.0 -1.1 -0.7  -10   15 1.6]]; % cWS Expanded 2010-11-27 CLZ; Moved left and up 2016-11-11 CLZ
B = [B;[ -5    4.0  6.0 -5.0 -3.6 -3.0  3.0 -1.1 -0.7  -35    5 1.6]]; % cSW Instances found, limits adjusted 2016-11-11 CLZ
B = [B;[  6    2.9  4.6  6.8  7.8 -3.0  3.0  0.7  1.1   80  120 1.6]]; % tWS CLZ 2007-12-12; Moved down, angle expanded left 2016-11-11 CLZ
B = [B;[ -6    5.2  7.7 -3.0 -0.7 -3.0  3.0  0.7  1.1  235  260 1.6]]; % tSW Expanded down, left 2016-11-11 CLZ
% = [B;[  8   -3.9 -1.5  7.7  9.5 -3.0  3.0  0.7  1.1  170  210 1.6]]; % tHh Expanded 2010-12-06 CLZ
B = [B;[  8   -3.9 -1.5  8.0 10.0 -3.0  3.0  0.7  1.1  150  200 1.6]]; % tHh Expanded 2010-12-06 CLZ; Moved up 2016-11-11 CLZ
B = [B;[  8.5  2.7  4.5  8.5  9.8 -3.0  3.0  0.7  1.1  175  200 1.6]]; % tHh New subcategory, different side of amino on C 2016-11-11 CLZ
B = [B;[  9   -5.3 -3.9  4.5  6.4 -3.0  3.0  0.7  1.1  -15   20 1.6]]; % cHS Tightened 2010-12-01 CLZ; Moved up and right 2016-11-11 CLZ
% = [B;[  9.1 -4.3 -1.0  6.5  8.8 -3.0  3.0  0.7  1.1  -70  -35 5.0]]; % cHS Shifted 2010-12-01 CLZ
B = [B;[  9.1 -4.3 -1.5  7.0  8.5 -3.0  3.0  0.7  1.1   25   65 1.6]]; % cHS Shifted 2010-12-01 CLZ; Moved down, tightened right 2016-11-11 CLZ
B = [B;[  9.2 -4.3 -3.5  6.2  7.0 -3.0  3.0  0.7  1.1   25   45 1.2]]; % cHS New, for bottom left corner of 9.1 2016-11-11 CLZ
% = [B;[ -9    3.2  5.3 -5.9 -4.5 -3.0  3.0  0.7  1.1  -40    5 1.6]]; % cSH
B = [B;[ -9    4.3  6.6 -5.2 -3.1 -3.0  3.0  0.7  1.1    5   25 1.6]]; % cSH Moved up and right 2016-11-11 CLZ
% = [B;[ 10   -3.4 -0.8  6.9  8.6 -3.0  3.0 -1.1 -0.7   70  105 1.6]]; % tHS
B = [B;[ 10   -2.5 -0.3  7.3  9.1 -3.0  3.0 -1.1 -0.7   65  105 1.6]]; % tHS Moved up and right, angle tightened 2016-11-11 CLZ
B = [B;[ 10.5 -5.1 -3.5  5.5  7.3 -3.0  3.0 -1.1 -0.7   90  120 1.6]]; % tHS New, more complete H-S overlap 2016-11-11 CLZ
B = [B;[ 11    2.7  4.6 -7.2 -5.0 -3.0  3.0 -1.1 -0.7  250  -75 2.6]]; % cSs Moved up, right, angle tightened left 2016-11-11 CLZ
% = [B;[-11    5.1  6.9 -3.7 -0.9 -3.3  3.3 -1.1 -0.7  -90  -30 1.6]]; % csS
B = [B;[-11    5.2  6.9 -3.8 -1.0 -3.3  3.3 -1.1 -0.7  -90  -60 1.6]]; % csS Tightened left, moved down, angle tightened right 2016-11-11 CLZ

s = size(B);
ClassLimits(1:s(1),1:s(2),14) = B;                      % CU is paircode 14

% GU pairs (paircode 15) ------------------------------------------------------

B =    [  1    5.0  6.9  6.0  7.5 -3.0  3.0 -1.1 -0.7   45   90 1.6];  % cWW, x widened, angle widened 2015-05-29 CLZ
B = [B;[  1.1  3.9  5.0  6.7  7.8 -3.0  3.0 -1.1 -0.7   60  100 1.6]]; % cWW
B = [B;[  1.2  5.5  7.5  4.7  6.0 -3.0  3.0 -1.1 -0.7   45   85 1.6]]; % cWW
% = [B;[  2    2.8  3.9  7.1  9.1 -3.4  3.0  0.7  1.1  180  215 5.0]]; % tWW
B = [B;[  2    2.6  5.2  7.1  9.0 -3.4  3.0  0.7  1.1  155  185 1.6]]; % tWW, x extended right and left, y top reduced, angle tightened 2015-05-29 CLZ
% = [B;[ -3   -6.1 -5.1  2.7  3.4 -3.0  3.0  0.7  1.1  -65  -35 5.0]]; % cHW
B = [B;[ -3   -6.1 -4.6  2.7  4.4 -3.0  3.0  0.7  1.1   40   65 1.6]]; % cHW, x extended right, y extended up 2015-05-29 CLZ
B = [B;[  4    5.1  7.9  6.5  8.2 -3.0  3.0 -1.1 -0.7  145  180 1.6]]; % tWH
% = [B;[ -4   -5.8 -4.8  5.1  6.1 -3.0  3.0 -1.1 -0.7  180  210 5.0]]; % tHW, x shifted right, angle extended left 2015-05-29 CLZ
B = [B;[ -4   -5.6 -3.8  5.1  6.1 -3.0  3.0 -1.1 -0.7  180  210 1.6]]; % tHW, x shifted right, angle extended left 2015-05-29 CLZ
% = [B;[  5    6.1  7.8  4.0  5.9 -3.0  3.0 -1.1 -0.7   15   40 5.0]]; % cWS
B = [B;[  5    6.6  7.9  4.0  5.9 -3.0  3.0 -1.1 -0.7    5   45 1.6]]; % cWS, x shifted to the right, ang extended left, right 2015-05-29 CLZ
B = [B;[  5.1  7.0  8.5  3.2  4.7 -3.0  3.0 -1.1 -0.7   30   55 1.6]]; % cWS
B = [B;[  5.2  7.7  8.7  3.7  4.8 -3.0  3.0 -1.1 -0.7   65   95 1.6]]; % cWS
% = [B;[ -5    2.8  3.8 -4.8 -3.8 -3.1  3.0 -1.1 -0.7  -65  -35 5.0]]; % cSW
B = [B;[ -5    2.8  3.9 -5.0 -4.0 -3.1  3.0 -1.1 -0.7  -60  -20 1.6]]; % cSW, ang shifted right, y shifted down, x extended right 2015-05-29 CLZ
% = [B;[  6    4.9  7.5  5.1  7.6 -3.0  3.0  0.7  1.1  235  -80 5.0]]; % tWS
B = [B;[  6    4.5  7.1  5.1  7.6 -3.0  3.0  0.7  1.1   80  130 1.6]]; % tWS, x shifted left, angle widened right 2015-05-30 CLZ
% = [B;[ -6    6.6  8.1 -3.7 -1.7 -3.0  3.0  0.7  1.1  115  155 5.0]]; % tSW
B = [B;[ -6    6.6  7.8 -3.4 -1.9 -3.0  3.0  0.7  1.1  205  240 1.6]]; % tSW, narrowed considerably to eliminate "pair" with stacked A in between 2015-06-03 CLZ
B = [B;[ -9    2.8  5.8 -6.2 -4.3 -3.0  3.0  0.7  1.1    0   30 1.7]]; % cSH
B = [B;[ -9.1  8.0  9.2 -4.7 -3.7 -3.0  3.0  0.7  1.1  -45  -15 1.6]]; % cSH
B = [B;[ -9.2  5.6  6.4 -5.3 -4.3 -3.0  3.0  0.7  1.1  -15    5 1.6]]; % cSH
B = [B;[-10    8.4 10.1 -3.7 -1.0 -3.0  3.0 -1.1 -0.7   50  100 1.6]]; % tSH
% = [B;[ 11    1.6  4.1 -8.7 -5.9 -3.0  3.0 -1.1 -0.7  230  -70 5.0]]; % cSs
B = [B;[ 11    1.6  4.1 -8.7 -7.2 -3.0  3.0 -1.1 -0.7  240  270 1.8]]; % cSs, vertical reduced, ang tightened right to eliminate stacked version of basepair 2015-06-03 CLZ
B = [B;[-11    5.7  7.5 -3.6 -1.9 -3.0  3.0 -1.1 -0.7  255  -30 1.6]]; % csS
B = [B;[-11.1  4.3  5.2 -4.0 -3.0 -3.0  3.0 -1.1 -0.7  250  270 1.6]]; % csS
% = [B;[ 12    8.1  8.9 -1.7  0.2 -3.0  3.0  0.7  1.1  180  220 5.0]]; % tSs
B = [B;[ 12    7.9  8.8 -1.8 -0.2 -3.0  3.0  0.7  1.1  150  180 1.6]]; % tSs, shifted down and left, ang moved left 2015-06-08 CLZ
% = [B;[ 13    7.5  9.2  5.6  7.5 -3.0  3.0 -1.1 -0.7  110  145 5.0]]; % bif
B = [B;[ 13    7.7  9.5  4.7  7.3 -3.0  3.0 -1.1 -0.7  105  140 1.6]]; % bif, extended down and right, ang move left 2015-06-08 CLZ

B = [B;[ 25   -5.4 -3.7  1.6  3.3 -3.0  3.0 -1.1 -0.7  125  165 1.6]]; % motif

s = size(B);
ClassLimits(1:s(1),1:s(2),15) = B;

% UU pairs (paircode 16) -----------------------------------------------------

B =    [  1    3.4  5.7  4.1  6.8 -3.0  3.0 -1.1 -0.7   35   80 1.6];  % cWw
% = [B;[  2    3.0  4.2  6.9  8.1 -3.4  3.0  0.7  1.1  170  200 5.0]]; % tWW 2013-12-04 moved from 2.1
B = [B;[  2    3.0  4.2  6.7  8.1 -3.4  3.0  0.7  1.1  165  195 1.6]]; % tWW 2013-12-04 moved from 2.1; vertical stretched down 2015-06-10 CLZ
% = [B;[  2.5  5.2  6.5  2.0  4.0 -3.0  3.0  0.7  1.1  160  200 5.0]]; % tWW 2008-10-23; moved from 2.1 to 2.5 2014-11-18
B = [B;[  2.5  5.2  6.4  2.0  3.9 -3.0  3.0  0.7  1.1  160  200 1.6]]; % tWW 2008-10-23; moved from 2.1 to 2.5 2014-11-18; upper right corner pulled in 2015-06-10 CLZ
% = [B;[  3    6.6  7.6  0.8  3.2 -3.0  3.0  0.7  1.1   55  100 5.0]]; % cWH Jesse_7_4_08
B = [B;[  3    6.5  8.2  0.7  3.9 -3.0  3.0  0.7  1.1  255  -55 1.6]]; % cWH Jesse_7_4_08; horizontal expanded, vertical expanded, ang shifted right 2015-06-10 CLZ
% = [B;[  4    3.8  6.2  5.2  7.1 -3.0  3.0 -1.1 -0.7  135  180 5.0]]; % tWH Expanded 2010-12-04 CLZ
B = [B;[  4    3.8  6.2  5.8  8.0 -3.0  3.0 -1.1 -0.7  140  175 1.6]]; % tWH Expanded 2010-12-04 CLZ; shifted up, ang tightened 2015-06-10 CLZ
% = [B;[  5    5.0  6.5  1.2  5.0 -3.0  3.0 -1.1 -0.7  -80   45 5.0]]; % cWS Jesse_4_18_08
B = [B;[  5    5.5  6.7  1.6  4.2 -3.0  3.0 -1.1 -0.7  -45   30 1.6]]; % cWS Jesse_4_18_08; shifted up and left, ang tightened left 2015-06-10 CLZ
% = [B;[  6    5.6  7.0  2.8  4.2 -3.0  3.0  0.7  1.1  200  230 5.0]]; % tWS
B = [B;[  6    3.8  6.6  3.3  5.9 -3.0  3.0  0.7  1.1   85  165 1.6]]; % tWS Expanded up and left, angle expanded right 2015-06-13 CLZ
B = [B;[  9   -6.6 -4.0  3.2  6.0 -3.0  3.0  0.7  1.1  -40   35 1.8]]; % cHS Expanded 2010-11-27 CLZ; a diffuse group, can't help it much 2015-06-13 CLZ
B = [B;[ 11    3.0  5.5 -6.5 -4.6 -3.0  3.0 -1.1 -0.7  260  -50 1.6]]; % cSs

s = size(B);
ClassLimits(1:s(1),1:s(2),16) = B;                      % UU is paircode 16

