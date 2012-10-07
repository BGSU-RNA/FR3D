% zDistanceToExemplars(Exemplar,Pair) computes the distance to each
% exemplar for the given pair of nucleotides

function [c,d,f,g,i] = zDistanceToExemplars(Exemplar,NT1,NT2)

pc = 4*(NT2.Code-1) + NT1.Code;     % AA is 1, CA is 2, etc.

k = 1;

F2.NT(1) = NT1;
F2.NT(2) = NT2;

for j = 1:length(Exemplar(:,pc)),
  E = Exemplar(j,pc);
  if ~isempty(Exemplar(j,pc).Class),
    c(k) = E.Class;
    d(k) = sqrt(xDiscrepancyFast(E,[NT1 NT2]))/2;
    f(k) = j;
    g(k) = pc;
  else
    c(k) = 99;
    d(k) = 99999999;
    f(k) = 1;               % fictitious
    g(k) = 1;               % fictitious
  end
  k = k + 1;
end

[a,i] = sort(d);

c = c(i);
d = d(i);
f = f(i);
g = g(i);
