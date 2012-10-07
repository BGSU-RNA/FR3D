% This function write rudimentary colors at the end of a PDB file

% Numbers is an nx1 vector telling how many nucleotides get each color
% Colors is a nx3 matrix of color specifiers

function [void] = nWriteColorInformation(fid,Numbers,Colors)

C = [];

for a = 1:length(Numbers),
  C = [C; ones(Numbers(a),1)*Colors(a,:)];
end

  % Write color information

  fprintf(fid,'SPDBVi 1 1 1 0 1 0 1 1 0 1 0 1 1 0 0\n');

  L = sum(Numbers);
  while L > 0,
    if L >= 20,
      fprintf(fid,'SPDBVf 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19\n');
      L = L - 20;
    else
      fprintf(fid,'SPDBVf ');
      for i = 1:L,
        fprintf(fid,'19 ');
      end
      fprintf(fid,'\n');
      L = 0;
    end
  end    

  % Write colors for each element

  L = sum(Numbers);
  a = 1;

  while L > 0,
    if L >= 4,
      fprintf(fid, 'SPDBVc ');
      fprintf(fid, '%4.2f ', C(a,:));
      fprintf(fid, '%4.2f ', C(a+1,:));
      fprintf(fid, '%4.2f ', C(a+2,:));
      fprintf(fid, '%4.2f ', C(a+3,:));
      fprintf(fid,'\n');
      L = L - 4;
      a = a + 4;
    else
      fprintf(fid, 'SPDBVc ');
      for i = 1:L,
        fprintf(fid,'%4.2f ', C(a,:));
        L = L - 1;
        a = a + 1;
      end
    end
  end

fprintf(fid,'\n');
