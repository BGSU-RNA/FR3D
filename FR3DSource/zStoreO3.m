% zStoreO3 adds the O3 prime atom of the previous nucleotide to the sugar of the
% current nucleotide, as row 13 for RNA sugars or as

function [File] = zStoreO3(File)

for f = 1:length(File),
	File(f).NT(1).Sugar(13,:) = File(f).NT(1).Sugar(10,:); % use phosphorus location
	if File(f).NT(1).Code == 9
		if isKey(File(f).NT(1).SugarDict, 'P')
			File(f).NT(1).SugarDict('O3''') = File(f).NT(1).SugarDict('P');
		end
	end
  for i = 2:length(File(f).NT),
    if strcmp(File(f).NT(i).Chain, File(f).NT(i-1).Chain), % same chain
	  PreviousO3 = File(f).NT(i-1).Sugar(5,:);
      	% PreviousO3 = File(f).NT(i-1).Sugar('O3''');

	  File(f).NT(i).Sugar(13,:) = PreviousO3;
      if File(f).NT(i).Code == 9
      	File(f).NT(i).SugarDict('PreviousO3''') = PreviousO3;
      end
    else
	  File(f).NT(i).Sugar(13,:) = [Inf Inf Inf]; % no previous O3'
	  if File(f).NT(i).Code == 9
		File(f).NT(i).SugarDict('PreviousO3''') = [Inf Inf Inf]; % no previous O3'
	  end
    end
  end
end
