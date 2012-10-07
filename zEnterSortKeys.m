% zEnterSortKeys(ViewParam) is part of zPairViewer.  It prompts the user
% for which fields to sort on, and in which order

function [ViewParam] = zEnterSortKeys(ViewParam)

fprintf('\n');
fprintf(' 1 - File number\n');
fprintf(' 2 - Paircode\n');
fprintf(' 3 - Displacement(1)\n');
fprintf(' 4 - Displacement(2)\n');
fprintf(' 5 - Displacement(3)\n');
fprintf(' 6 - Normal(3)\n');
fprintf(' 7 - Angle between planes\n');
fprintf(' 8 - Angle of rotation\n');
fprintf(' 9 - Absolute value of Gap\n');
fprintf('10 - Computer classification\n');
fprintf('11 - Hand classification\n');
fprintf('12 - Minimum distance\n');
fprintf('13 - Index of first nucleotide\n');
fprintf('14 - Lowest index in pair\n');
fprintf('15 - Angle in the plane of the first base\n')
fprintf('16 - Smallest hydrogen bond angle\n')
fprintf('17 - Largest hydrogen bond angle\n')
fprintf('18 - Smallest hydrogen bond length\n')
fprintf('19 - Largest hydrogen bond length\n')
fprintf('20 - Distance to nearest exemplar\n')
fprintf('21 - Classification of nearest exemplar\n')
fprintf('22 - Degree of stacking overlap\n')
fprintf('23 - Discrepancy with specified pair\n')

fprintf('\n');

fprintf('Enter sort key(s) in brackets, negative for reversed order\n');
inp = input('');
ViewParam.SortKeys = inp; 
