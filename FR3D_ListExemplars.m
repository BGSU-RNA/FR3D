% zListExamplars lists all basepair exemplars

pcodes = [1 5 6 7 9 11 13 14 15 16];
% 1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

% load exemplars from previous session -------------------------------------

 load('PairExemplars','Exemplar');

% loop through paircodes and computer classifications ----------------------

Entries = [];
c = 1;

for j = 1:length(pcodes),
 Paircode = pcodes(j);
 for row = 1:length(Exemplar(:,Paircode)),

  E = Exemplar(row,Paircode);

  if ~isempty(E.Filename),
    n1 = E.NT1;
    n2 = E.NT2;
    
    a = [sprintf('%s%s ',  n1.Base, n2.Base) ... 
         sprintf('%s', E.Pair.EdgeText) ...
         sprintf('%10s %2s %5s ',  E.Filename, n1.Base, n1.Number) ... 
         sprintf('%2s %5s ', n2.Base, n2.Number) ...
         sprintf('%6.2f ', E.Class) ...
         sprintf('%5d instances\n', E.Count)];

    Entries{c} = a;
    b(c) = E.Class;
    c = c + 1;

  end
 end
end

[y,i] = sort(abs(b));
for j = 1:length(Entries),
  fprintf('%s', Entries{i(j)});
end

