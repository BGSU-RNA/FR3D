% zWriteClassLimitsPython.m reads zClassLimits.m and writes out the limits as Python code
% It simply uses the geometry of RNA basepairs and applies them to DNA basepairs
% We'll see how that works

[ClassLimits,CurrentVersion] = zClassLimits;

[NumCategories,Columns,Codes] = size(ClassLimits);

Paircode{ 1} = 'DA,DA';
Paircode{ 5} = 'DA,DC';
Paircode{ 6} = 'DC,DC';
Paircode{ 7} = 'DG,DC';
Paircode{ 9} = 'DA,DG';
Paircode{11} = 'DG,DG';
Paircode{13} = 'DA,DT';
Paircode{14} = 'DC,DT';
Paircode{15} = 'DG,DT';
Paircode{16} = 'DT,DT';

Basecode{ 1} = [1 1];
Basecode{ 5} = [1 2];
Basecode{ 6} = [2 2];
Basecode{ 7} = [3 2];
Basecode{ 9} = [1 3];
Basecode{11} = [3 3];
Basecode{13} = [1 4];
Basecode{14} = [2 4];
Basecode{15} = [3 4];
Basecode{16} = [4 4];

ot = 'nt_nt_cutoffs = {}\n';
ot = '';                        % this can be run immediately after RNA limits

for k = 1:16
	pc = Paircode{k};
	if length(pc) > 0
		SubcategoryCounter = 0;
		PreviousCategory = 0;
		PreviousCategoryText = '';
		ot = [ot sprintf('nt_nt_cutoffs["%s"] = {}\n',pc)];

		for j = 1:NumCategories
			if abs(ClassLimits(j,1,k)) > 0 && abs(ClassLimits(j,1,k)) < 14
				CL = ClassLimits(j,:,k);
				CurrentCategory = fix(CL(1));
				CurrentSubCategory = abs(CL(1) - CurrentCategory);
				CategoryText = strtrim(zEdgeText(CurrentCategory,0,Basecode{k}(1),Basecode{k}(2)));
				if CurrentCategory == 13
					CategoryText = 'cWB';
				end
				if CurrentCategory == -13
					CategoryText = 'cBW';
				end

				if CurrentSubCategory >= 0.5
					CategoryText = "a" + CategoryText;      % "alternative", let's see how it works
				end

				if CurrentCategory ~= PreviousCategory || ~strcmp(CategoryText,PreviousCategoryText)
					SubcategoryCounter = 0;
					ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"] = {}\n',pc,CategoryText)];
				end

				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d] = {}\n',pc,CategoryText,SubcategoryCounter)];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["xmin"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(2))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["xmax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(3))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["ymin"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(4))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["ymax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(5))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["zmin"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(6))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["zmax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(7))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["normalmin"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(8))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["normalmax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(9))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["anglemin"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(10))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["anglemax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(11))];
				ot = [ot sprintf('nt_nt_cutoffs["%s"]["%s"][%d]["gapmax"] = %0.2f\n',pc,CategoryText,SubcategoryCounter,CL(12))];

				SubcategoryCounter = SubcategoryCounter + 1;
				PreviousCategory = CurrentCategory;
				PreviousCategoryText = CategoryText;

			end
		end
	end
end

fprintf(ot)
