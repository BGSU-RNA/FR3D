% xExcludeOverlap(Candidates) removes candidates which overlap with one
% another, leaving only the one that occurs first in the list

function [Candidates,Discrepancy] = xExcludeOverlap(Candidates,Discrepancy,Limit)

N      = length(Candidates(1,:)) - 1;  % number of nucleotides

OK     = zeros(size(Candidates(:,1))); % rows of Candidates which are kept
OK(1)  = 1;                            % keep the first; lowest discrepancy
NumOK  = 1;                            % counter for number kept

for i = 2:length(Candidates(:,1)),     % 
  AddCand = 1;                         % default is to add this one
  if Discrepancy(i) >= 0,           % if discrepancy is below GuarCutoff
    j = 1;                               % first kept candidate
    while (j <= NumOK) & (AddCand == 1), % go through kept candidates
      if Candidates(i,N+1) == Candidates(OK(j),N+1),  % same file
        Both = intersect(Candidates(i,1:N), Candidates(OK(j),1:N)); %overlap
        if length(Both) > N/2,         % lots of overlap here
          AddCand = 0;                 % end the loop
        end
      end
      j = j + 1;                       % on to next kept candidate
    end
    if AddCand == 1,                   % no overlapping candidate found
      NumOK = NumOK + 1;
      OK(NumOK) = i;
    end
  end

  if NumOK >= Limit,
    fprintf('Top %2d distinct candidates retained\n', Limit);
    break
  end

end

OK = OK(1:NumOK);

Candidates = Candidates(OK,:);
Discrepancy = Discrepancy(OK);
