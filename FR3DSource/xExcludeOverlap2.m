% xExcludeOverlap2(Candidates,Discrepancy,Perm) removes candidates which overlap with one another, leaving only the one that occurs first in the list

function [Candidates,Discrepancy] = xExcludeOverlap2(Candidates, Discrepancy, Perm)

N    = length(Candidates(1,:)) - 1;     % number of nucleotides

Keep = ones(size(Candidates(:,1)));     % rows of Candidates which are kept

[C,ord] = sortrows(Candidates,[N+1 Perm]); % sort by filenumber, nucleotides

%D = diff(C);                           % matrix of row differences
%D(:,1:(N+1)) = abs(sign(D(:,1:(N+1))));% 1 for every difference, 0 if same

a = 1;                                  % index of current lowest discrepancy
for i = 2:(length(C(:,1))-1),           % go through row differences
  if C(a,N+1) == C(i,N+1),              % if file numbers agree
    if sum(C(a,1:N) == C(i,1:N)) > N/2, % if more than half agree
        if Discrepancy(ord(i)) > Discrepancy(ord(a)),
        Keep(i) = 0;                    % reject i; a has lower discrep
      else
        Keep(a) = 0;                    % reject a; i has lower discrep
        a = i;                          % found one with lower discrepancy
      end
    else
      a = i;                            % moved on to new motif
    end
  else
    a = i;                              % moved on to new file
  end
end

Keep(ord) = Keep;                       % undo the re-ordering
OK = find(Keep);
Candidates = Candidates(OK,:);         % return only non-overlapping candidates
Discrepancy = Discrepancy(OK);
