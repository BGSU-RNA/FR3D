% zExemplarTablesExcel collects tables of discrepancies formatted for Excel

TT = {};
UU = {};

for c = 1:12,
  [T,U] = zExemplarTable(c,1,0,0);

  [a,b] = size(T);
  [c,d] = size(TT);
  TT((c+2):(c+1+a),1:b) = T;

  [a,b] = size(U);
  [c,d] = size(UU);
  UU((c+2):(c+1+a),1:b) = U;
end

Lett = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% ---------------------------------------------- Write 16x16 tables

xlswrite('c:\IsoDiscrepancyTables2.xls',TT);
Excel = actxserver('Excel.Application');
Excel.Workbooks.Open('c:\IsoDiscrepancyTables2.xls');

[c,d] = size(TT);
for i = 1:c,
  Range = Excel.Range(['AA' num2str(i)]);
  Range.Interior.ColorIndex = min(56,i);
  for j = 6:d,
    if strcmp(class(TT{i,j}),'double') && ~isempty(TT{i,j}),
      Range = Excel.Range([Lett(j) num2str(i)]);
      if TT{i,j} < 2,
        Range.Interior.ColorIndex = 3;
      elseif TT{i,j} < 3.3,
        Range.Interior.ColorIndex = 6;
      elseif TT{i,j} < 5.5,
        Range.Interior.ColorIndex = 8;
      else 
        Range.Interior.ColorIndex = 5;
      end
    end
  end
end

%Range.Font.ColorIndex = 42;
Excel.Visible=1;

% ---------------------------------------------- Write 4x4 tables

xlswrite('c:\IsoDiscrepancyTables3.xls',UU);
Excel = actxserver('Excel.Application');
Excel.Workbooks.Open('c:\IsoDiscrepancyTables3.xls');

[c,d] = size(UU);
for i = 1:c,
  for j = 2:d,
    if strcmp(class(UU{i,j}),'double') && ~isempty(UU{i,j}),
      Range = Excel.Range([Lett(j) num2str(i)]);
      if UU{i,j} < 2,
        Range.Interior.ColorIndex = 3;
      elseif UU{i,j} < 3.3,
        Range.Interior.ColorIndex = 6;
      elseif UU{i,j} < 5.5,
        Range.Interior.ColorIndex = 8;
      else 
        Range.Interior.ColorIndex = 5;
      end
    end
  end
end

%Range.Font.ColorIndex = 42;
Excel.Visible=1;
