% zDisplayPairCounts displays counts of pairwise interactions, either from exemplars or from structures.  The format is Count(base1code,base2code,category).

function [void] = zDisplayPairCounts(Count,CountFilename)

Total = [];

[s,t,MaxCat] = size(Count);

if MaxCat > 100,
  m = 10;                           % category numbers are all multiplied by 10
else
  m = 1;
end

for ca = 1:12,
  Total(ca) = sum(sum(Count(:,:,m*ca)));
end

for ca = [1 2 7 8];
  for c1 = 2:4,
    for c2 = 1:(c1-1);
%      Count(c1,c2,m*ca) = -Count(c2,c1,m*ca);
    end
  end
end

row = 1;
clear T

T{row,1} = 'Family';
T{row,2} = 'Count';
T{row,3} = '% of total';

for ca = 1:12,                           % families 1-12 only
  T{row+ca,1} = zEdgeText(ca,0);
  T{row+ca,2} = Total(ca);
  T{row+ca,3} = 100*Total(ca)/sum(Total(1:12));
end

row = row + ca + 2;

Letters = 'ACGU';

S = 6;

for ca = 1:12,
  T{row,1} = zEdgeText(ca,0);
  T{row,1+S} = zEdgeText(ca,0);
  for i = 1:4,
    T{row+i,1} = Letters(i);
    T{row,1+i} = Letters(i);
    for j = 1:4,
      T{row+j,1+i} = Count(j,i,m*ca);
    end

    T{row+i,1+S} = Letters(i);
    T{row,1+i+S} = Letters(i);
    for j = 1:4,
      T{row+j,1+i+S} = 100*Count(j,i,m*ca)/Total(ca);
    end

  end

  row = row + 6;
end

T{row,2} = 'Counts from 3D structures';
row = row + 1;
T{row,18} = 'Total';
T{row,19} = '% of total';
for i = 1:4,
  for j = 1:4,
    pc = 4*(j-1) + i;
    T{row,1+pc} = [Letters(j) Letters(i)];
  end
end

for ca = 1:12,
  T{row+ca,1} = zEdgeText(ca,0);
  for i = 1:4,
    for j = 1:4,
      pc = 4*(j-1) + i;
      if Count(i,j,m*ca) ~= 0,
        T{row+ca,1+pc} = Count(j,i,m*ca);
      end
    end
  end
  T{row+ca,18} = Total(ca);
  T{row+ca,19} = 100*Total(ca)/sum(Total(1:12));
end

row = row + 15;

T{row,2} = 'Frequencies from 3D structures';
row = row + 1;
T{row,18} = 'Total';
T{row,19} = '% of total';
for i = 1:4,
  for j = 1:4,
    pc = 4*(j-1) + i;
    T{row,1+pc} = [Letters(j) Letters(i)];
  end
end

for ca = 1:12,
  T{row+ca,1} = zEdgeText(ca,0);
  for i = 1:4,
    for j = 1:4,
      pc = 4*(j-1) + i;
      if Count(i,j,m*ca) ~= 0,
        T{row+ca,1+pc} = 100*Count(j,i,m*ca)/Total(ca);
      end
    end
  end
  T{row+ca,18} = 100;
  T{row+ca,19} = 100*Total(ca)/sum(Total(1:12));
end

xlswrite(CountFilename,T);

