% xMutualDiscrepancy computes the mutual distances between candidates

function [Search] = xMutualDiscrepancy(File,Search)

Query      = Search.Query;
Candidates = Search.Candidates;

N     = Query.NumNT;
[s,t] = size(Candidates);

if ~isfield(Query,'LocWeight'),              % should not be necessary
  Query.LocWeight = ones(1,N);
end  

if ~isfield(Query,'AngleWeight'),            % should not be necessary
  Query.AngleWeight = ones(1,N);
end 

% ----------------------------------- Compute mutual discrepancies btw cand's

if min(Search.DiscComputed(1,1:s)) == 0,             % some not computed
  for k=1:s,
    if Search.DiscComputed(k) == 0,                % k has never been done
      f1             = Candidates(k,N+1);
      c1.NumNT       = Query.NumNT;                  % set up c1 as "Model"
      c1.NT          = File(f1).NT(Candidates(k,1:N));
      c1.LocWeight   = Query.LocWeight;
      c1.AngleWeight = Query.AngleWeight;
      c1.LDiscCutoff = Inf;
      c1 = xPrecomputeForDiscrepancy(c1);

      for j=1:k-1,
        f2 = Candidates(j,N+1);
        c2 = File(f2).NT(Candidates(j,1:N));     
        Search.Disc(j,k) = sqrt(xDiscrepancyFast(c1,c2))/Query.NumNT;
        Search.Disc(k,j) = Search.Disc(j,k);
      end
      Search.DiscComputed(k) = 1;
    end
  end
end

