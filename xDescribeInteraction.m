
clear Message

Message{1} = 'Basepair, base stacking, or letter pair constraint is specified above the diagonal.';

Message{2} = 'To specify that all candidate motifs must have a tWH basepair between the nucleotides corresponding to the first and second nucleotides in the query motif, type tWH in the first row, second column.  This means that the nucleotide in the first row must use its Watson-Crick edge, and the nucleotide in the second column must use its Hoogsteen edge.';

Message{3} = 'Valid basepair specifications are cWW, tWW, cWH, cHW, tWH, tHW, cWS, cSW, tWS, tSW, cHH, tHH, cHS, cSH, tHS, tSH, cSS, tSS. Note, however, that the cSS and tSS interactions are not, in fact, symmetric, because each base can use the sugar edge differently. Follwing Leontis, Stombaugh, Westhof (NAR 2002), type cSs to specify that the first base has priority, csS for the second, or cSS for either.';

Message{4} = 'Specifying multiple interactions allows more ways a candidate can satisfy the constraints; for example, typing cWH cHW requires a cis Watson-Crick/Hoogsteen basepair, but either base can use the Watson-Crick edge, and the other uses the Hoogsteen edge. ';

Message{5} = 'The abbreviation "trans" gives all trans categories, "cis" for cis. ';

Message{6} = 'Type "bif" for bifurcated basepairs (see LSW 2002). ';

Message{7} = 'Type ~cWW to exclude candidates having a cWW basepair. ';

Message{8} = 'Some pairs of bases are close to, say, cWW, but do not meet the strict criteria for membership in the cWW classification. Type ncWW ("near cWW") to get basepairs that are  not classified into any category, but for which the cWW category is the closest match, up to a certain fairly generous limit. Type "cWW ncWW" to get cWW and near cWW pairs, cWW. Type ntrans to get all pairs nearest to a trans pair. ';

Message{9} = 'Type s35 for stacking in which the first base uses its 3 face, and the second base uses its 5 face. Similarly, type s53, s33, or s55. Type "stack" to allow all stacking interactions. The prefixes "n" and "~" work with stacking, as above. ';

Message{10} = 'To specify that the nucleotides must match a certain pattern, type, for example, "cWW CG GC" to get only CG or GC cWW pairs.';

mEditBoxWrap(Message,'Interaction help');

