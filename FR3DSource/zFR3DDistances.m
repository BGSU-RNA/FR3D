% zFR3DDistances calculate pairwise distances between nucleotides for a FR3D search
% The largest distance recorded is Query.DistCutoff, except when the Flank constraint (aka BorderSS) is active, in which case
% it records distances up to a larger limit to accommodate all single-stranded regions

function File = zFR3DDistances(File,Query)

BorderSSActive = 0;
[s,t] = size(Query.Flank);
for a = 1:s,
	for b = 1:t,
		if ~isempty(Query.Flank),
			BorderSSActive = 1;
		end
	end
end

for f = 1:length(File),

  c = cat(1,File(f).NT(1:File(f).NumNT).Center);

  if BorderSSActive,
  	[i,j] = find(triu(File(f).Flank));                     % pairs of nucleotides which satisfy the BorderSS relation
  	maxdist = 0;
  	for k = 1:length(i),
  		maxdist = max(maxdist,zDistance(c(i(k),:),c(j(k),:)));
  	end
	  File(f).Distance = zMutualDistance(c,max(Query.DistCutoff,maxdist+10)); 

  else
	  File(f).Distance = zMutualDistance(c,Query.DistCutoff);    % sparse matrix of center-center distances, up to Query.DistCutoff
  end
end
