% zDistanceToExemplar(Exemplar,NT1,NT2) computes the distance between the pair NT1, NT2 of exemplars that have been classified as c to the exemplar for that interaction

function [d] = zDistanceToExemplars(Exemplar,N1,N2,c)

paircode = 4*(N2.Code-1) + N1.Code;           % AA is 1, CA is 2, etc.

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

% change order of nucleotides, if needed

switch paircode
  case {2, 3, 4, 8, 10, 12},                  % put N2 at the origin
    M1 = N2;
    M2 = N1;
    s  = -1;                                  % bases in reversed order
    if paircode ~= 10,
      c  = -c;
    end
    pc = 4*(M2.Code-1) + M1.Code;             % calculate new paircode
  otherwise
    M1 = N1;
    M2 = N2;
    if paircode == 7,
      c = -c;
    end
    s  = 1;                                   % bases in original order
    pc = paircode;
end

e = find(cat(2,Exemplar(:,pc).Class) == c);

if ~isempty(e),
  E = Exemplar(e(1),pc);
  d = sqrt(xDiscrepancyFast(E,[M1 M2]))/2;
  if M1.Code == M2.Code,
    d = min(d,sqrt(xDiscrepancyFast(E,[M2 M1]))/2);
  end
elseif M1.Code == M2.Code,
  e = find(cat(2,Exemplar(:,pc).Class) == -c);
  E = Exemplar(e(1),pc);
  d = sqrt(xDiscrepancyFast(E,[M1 M2]))/2;
  if M1.Code == M2.Code,
    d = min(d,sqrt(xDiscrepancyFast(E,[M2 M1]))/2);
  end
else
  e = find(fix(cat(2,Exemplar(:,pc).Class)) == fix(c));
  if ~isempty(e),  
    E = Exemplar(e(1),pc);
    d = sqrt(xDiscrepancyFast(E,[M1 M2]))/2;
  else
    fprintf('Pair %s%4s - %s%4s %s distance Inf to exemplar\n', M1.Base,M1.Number,M2.Base,M2.Number, zEdgeText(c));
    d = Inf;
  end
end

