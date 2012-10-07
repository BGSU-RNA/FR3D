
clear Message

Message{1} = 'Nucleotide mask, location weight, and angle weight are specified  in the boxes on the diagonal.';

Message{2} = 'For the mask, type A to allow only A''s, similarly with C, G, U. Type AC to allow only A or C, or AGU to allow only A, G, or U. Type ~G to exclude G (or type ACU). You may use one of these standard abbreviations:   M = AC;  R = AG;  W = AU;  S = CG;  Y = CU;  K = GU   V = ACG; H = ACU; D = AGU; B = CGU; N = ACGU Using no mask is the same as specifying ACGU or N. Self-contradictory specifications (such as A ~A) may give unpredictable results or empty searches.';

Message{3} = 'Specify the location weight for the given nucleotide by typing something of the form LW2.5 in the box. Location weights are normalized automatically, so making them all equal to 5.0 has no effect on the search. The default value, if nothing is specified, is 1.0.';

Message{4} = 'Specify the angle weight for the given nucleotide by typing something of the form AW2.0 in the box. Angle weights are not automatically normalized; the default is 1.0. Making angle weights larger than 1.0 increases the orientation error. Thus, candidates with low discrepancy will have bases more closely aligned than when using the default value of 1.0 for angle weight.';

Message{5} = 'Separate the mask, location weight, and angle weight by spaces. ';

mEditBoxWrap(Message,'Mask help');

