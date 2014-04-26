% zExemplarFrequencyTables loads the file of pair exemplars and makes tables of occurrence frequencies

% column numbers  1-AA  2-AC  3-AG  4-AU  5-CA  6-CC  7-CG  8-CU  9-GA 10-GC 11-GG 12-GU 13-UA 14-UC 15-UG 16-UU

function [void] = zExemplarFrequencyTables(T)
	
RNA = 'ACGU';

for c1 = 1:4,
	for c2 = 1:4,
		cn  = 4*(c1-1)+c2;          % column number
		bc = [RNA(c1) RNA(c2)];
		fprintf('%2d-%s ',cn,bc);
	end
end
fprintf('\n');

for t = 1:12,             % loop through basepair families
	for c1 = 1:4,
		for c2 = 1:4,
			cn  = 4*(c1-1)+c2;          % column number
			bc = [RNA(c1) RNA(c2)];
			j = find(ismember(upper(T{t}(:,2)),bc));
			if ~isempty(j),
				Count(t,cn) = str2num(T{t}{j(1),7});
				Found(t,cn) = 1;
			else
				Count(t,cn) = 0;
				Found(t,cn) = 0;
			end
		end
	end
end

% remove double counting in symmetric families

C = Count;
repeatcolumns = [5 9 10 13 14 15];
symmetricfamilies = [1 2 7 8];
C(symmetricfamilies,repeatcolumns) = 0;
nonrepeatcolumns = setdiff(1:16,repeatcolumns);

fid = fopen(['Exemplars' filesep 'Exemplars_12x1_table.txt'],'w');
for t = 1:12,
	bpf = zEdgeText(t);
	bpf = strrep(bpf,' ','');
	fprintf(fid,'%s\t%d\t%6.2f\n',bpf,sum(C(t,:)),100*sum(C(t,:))/sum(sum(C)));
end
fclose(fid);

fid = fopen(['Exemplars' filesep 'Exemplars_4x4_table.txt'],'w');
fprintf(fid,'\tA\tC\tG\tU\n');
for c1 = 1:4,
	fprintf(fid,'%s\t',RNA(c1));
	for c2 = 1:4,
		cn  = 4*(c1-1)+c2;          % column number
		if c2 > c1,
			cn = [cn 4*(c2-1)+c1];          % two column numbers
		end
		if c2 >= c1,
			fprintf(fid,'%6.2f',100*sum(sum(C(:,cn)))/sum(sum(C)));		
		end
		if c2 < 4,
			fprintf(fid,'\t');
		end
	end
	fprintf(fid,'\n');
end
fclose(fid);

fid = fopen(['Exemplars' filesep 'Exemplars_12x16_table.txt'],'w');
fprintf(fid,'\t');
for c1 = 1:4,
	for c2 = 1:4,
		bc = [RNA(c1) RNA(c2)];
		fprintf(fid,'%s',bc);
		if c1 < 4 || c2 < 4,
			fprintf(fid,'\t');
		else
			fprintf(fid,'\n');
		end
	end
end
for t = 1:12,
	bpf = zEdgeText(t);
	bpf = strrep(bpf,' ','');
	fprintf(fid,'%s\t',bpf);
	for cn = 1:16,
		if max(t == symmetricfamilies) == 0 || max(cn == repeatcolumns) == 0,
			if Found(t,cn) == 0,
				fprintf(fid,'NA');
			else
				fprintf(fid,'%6.2f',100*C(t,cn)/sum(C(t,:)));
			end
		else
			fprintf(fid,'*');
		end
		if cn < 16,
			fprintf(fid,'\t');
		else
			fprintf(fid,'\n');
		end
	end
end

fid = fopen(['Exemplars' filesep 'Exemplars_12x10_table.txt'],'w');
fprintf(fid,'\t');
for c1 = 1:4,
	for c2 = 1:4,
		cn  = 4*(c1-1)+c2;          % column number
		if max(cn == nonrepeatcolumns) == 1,
			bc = [RNA(c1) RNA(c2)];
			fprintf(fid,'%s',bc);
			if c1 < 4 || c2 < 4,
				fprintf(fid,'\t');
			else
				fprintf(fid,'\n');
			end
		end
	end
end
for t = 1:12,
	bpf = zEdgeText(t);
	bpf = strrep(bpf,' ','');
	fprintf(fid,'%s\t',bpf);
	for c1 = 1:4,
		for c2 = 1:4,
			cn  = 4*(c1-1)+c2;          % column number
			if max(cn == nonrepeatcolumns) == 1,
				if c1 ~= c2,
					cn = [cn 4*(c2-1)+c2];
				end
				fprintf(fid,'%6.2f',100*sum(C(t,cn))/sum(sum(C(:,cn))));
				if c1 < 4 || c2 < 4,
					fprintf(fid,'\t');
				else
					fprintf(fid,'\n');
				end
			end
		end
	end
end
fclose(fid);
