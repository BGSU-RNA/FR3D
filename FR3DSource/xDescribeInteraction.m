
clear Message

Message{1} = 'Basepair, base stacking, base phosphate, or letter pair constraint is specified above the diagonal.';

Message{2} = '';

Message{3} = 'To specify that all candidate motifs must have a tWH basepair between the nucleotides corresponding to the first and second nucleotides in the query motif, type tWH in the first row, second column.  This means that the nucleotide in the first row must use its Watson-Crick edge, and the nucleotide in the second column must use its Hoogsteen edge.';

Message{4} = '';

Message{5} = 'Valid basepair specifications are cWW, tWW, cWH, cHW, tWH, tHW, cWS, cSW, tWS, tSW, cHH, tHH, cHS, cSH, tHS, tSH, cSS, tSS. Note, however, that the cSS and tSS interactions are not, in fact, symmetric, because each base can use the sugar edge differently. Follwing Leontis, Stombaugh, Westhof (NAR 2002), type cSs to specify that the first base has priority, csS for the second, or cSS for either.';

Message{6} = '';

Message{7} = 'Specifying multiple interactions allows more ways a candidate can satisfy the constraints; for example, typing cWH cHW requires a cis Watson-Crick/Hoogsteen basepair, but either base can use the Watson-Crick edge, and the other uses the Hoogsteen edge. ';

Message{8} = '';

Message{9} = 'The abbreviation "trans" gives all trans categories, "cis" for cis. ';

Message{10} = '';

Message{11} = 'Type "bif" for bifurcated basepairs (see LSW 2002). ';

Message{12} = '';

Message{13} = 'Type ~cWW to exclude candidates having a cWW basepair. ';

Message{14} = '';

Message{15} = 'Some pairs of bases are close to, say, cWW, but do not meet the strict criteria for membership in the cWW classification. Type ncWW ("near cWW") to get basepairs that are  not classified into any category, but for which the cWW category is the closest match, up to a certain fairly generous limit. Type "cWW ncWW" to get cWW and near cWW pairs, cWW. Type ntrans to get all pairs nearest to a trans pair. ';

Message{16} = '';

Message{17} = 'Type s35 for stacking in which the first base uses its 3 face, and the second base uses its 5 face. Similarly, type s53, s33, or s55. Type "stack" to allow all stacking interactions. The prefixes "n" and "~" work with stacking, as above. ';

Message{18} = '';

Message{19} = 'To specify that the nucleotides must match a certain pattern, type, for example, "cWW CG GC" to get only CG or GC cWW pairs.';

Message{20} = '';

Message{21} = 'To require that two nucleotides make a base-phosphate interaction, enter BPh in the corresponding yellow box. This will select pairs of nucleotides in which the first nucleotide''s base is a hydrogen bond donor and the second nucleotide''s phosphate is an acceptor. To reverse the roles, type PhB. To specify particular base-phosphate categories, type 0BPh, 1BPh, 2BPh, ..., 9BPh, or 0PhB, 1PhB, etc. For near base-phosphate interactions, type nBPh, nPhB, n1BPh, n1PhB, etc.';

Message{22} = '';

Message{23} = 'One can restrict to pairs that play a certain role in the secondary and tertiary structure.  For pairs that are nested, type "N" or "nested".  For pairs that cross nested interactions but involve nucleotides in the same branch of the RNA, type "local" or "L".  For long-range or distant interactions, between different branches of the RNA, type "long-range", "distant", "D", or "LR".  Note that "nested", "local", and "distant" are mutually exclusive.  They can be negated with ~, but ~local only returns distant interactions, not nested ones.';

Message{24} = '';

Message{25} = 'To specify bases on the same strand that form cWW pairs that flank a hairpin, internal, or junction loop, type "flank" or "F".  Note:  for internal and junction loops, flanking nucleotides will be on the same strand, one on each side of the internal or junction loop.  Such flanking nucleotides usually do not interact with one another.';

mEditBoxWrap(Message,'Interaction help');

