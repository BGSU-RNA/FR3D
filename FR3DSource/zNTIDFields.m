% zNTIDFields breaks apart the text of a nucleotide ID into a cell array of fields, some of which may be blank

function [fields] = zNTIDFields(s)

fields{1} = '';
fields{2} = '';
fields{3} = '';
fields{4} = '';
fields{5} = '';
fields{6} = '';
fields{7} = '';
fields{8} = '';
fields{9} = '';

if length(s) > 0
	j = find(s == '|');
	fields{1} = s(1:(j(1)-1));
	for k = 1:(length(j)-1),
		fields{k+1} = s((j(k)+1):(j(k+1)-1));
	end
	fields{length(j)+1} = s((j(end)+1):end);

	if length(fields) < 9 || isempty(fields{9}),
		fields{9} = '1_555';
	end
end
