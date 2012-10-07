% zWriteHandFile writes Filename.hand

function [void] = zWriteHandFile(File)

[i, j] = find(File.CI > 0);                     % find pairs with hand class
k = find(i<j);                                  % look at each pair only once
i = i(k);                                       % reduce list of indices
j = j(k);                                       % reduce list of indices

if length(i) > 0,
  fid = fopen([File.Filename '.hand'],'w');       % open for writing
  for n=1:length(i),
    p = i(n);                                     % get first base index
    q = j(n);                                     % and second base index
    ci = File.CI(p,q);                                 % hand class index
    if (abs(File.HandClass(ci)) > 0) || (length(File.Comment{ci}) > 0),
      fprintf(fid,'%2s ', File.NT(p).Base);
      fprintf(fid,'%5s ', File.NT(p).Number);
      fprintf(fid,'%2s ', File.NT(p).Chain);
      fprintf(fid,'%2s ', File.NT(q).Base);
      fprintf(fid,'%5s ', File.NT(q).Number);
      fprintf(fid,'%2s ', File.NT(q).Chain);
      fprintf(fid,'%6.2f ', File.Inter(p,q));
      fprintf(fid,'%6.2f ', File.HandClass(ci));
      if length(File.Comment{ci}) > 0,
        fprintf(fid,'%s',File.Comment{ci});
      end
      fprintf(fid,'\n');
    end
  end

  fclose(fid);
  fprintf('Wrote %s.hand\n',File.Filename);
end
