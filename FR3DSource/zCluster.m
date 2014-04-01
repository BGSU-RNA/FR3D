% zCluster(D) returns a listing of groupings for hierarchical clustering of data points with mutual distances in matrix D

function [Z] = zCluster(D)

	[a,b] = size(D);
	Z = [];

	if a == b,
		Groups = {};
		for k = 1:a,
			D(k,k) = Inf;
			Groups{k} = k;
		end
		RowToGroup = 1:a;        % mapping from current row to group

		GD = D;                  % group distances, will shrink with each group merger

		while length(GD(1,:)) > 1,
			[y,c] = min(GD);
			[z,r] = min(y);
			r = r(1);              % these are the two groups to merge
			c = c(r);

      if 0 > 1,
      	GD
      	[r c z]
      end

			Z = [Z; [RowToGroup(r) RowToGroup(c) z]]; % Matlab's way of recording which groups were merged

			k = length(Groups) + 1;
			Groups{k} = [Groups{r} Groups{c}];    % merge the groups

			Remaining = setdiff(1:a,[r c]);       % remove these two

			newD = min(GD(Groups{k},Remaining));  % single linkage

			RowToGroup = [RowToGroup(Remaining) k];
			GD = [GD(Remaining,Remaining) newD'; newD Inf];

			[a,b] = size(GD);

		end

	else
		fprintf('zCluster cannot work with a non-square matrix\n');
	end