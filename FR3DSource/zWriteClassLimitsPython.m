% zWriteClassLimitsPython.m reads zClassLimits.m and writes out the limits as Python code

[ClassLimits,CurrentVersion] = zClassLimits;

[NumCategories,Columns,Codes] = size(ClassLimits);

Paircode{ 1} = 'A,A';
Paircode{ 5} = 'A,C';
Paircode{ 6} = 'C,C';
Paircode{ 7} = 'G,C';
Paircode{ 9} = 'A,G';
Paircode{11} = 'G,G';
Paircode{13} = 'A,U';
Paircode{14} = 'C,U';
Paircode{15} = 'G,U';
Paircode{16} = 'U,U';

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
