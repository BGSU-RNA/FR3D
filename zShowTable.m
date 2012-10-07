
function [void] = showtable(rowlabels, collabels, table, Title)

if nargin < 4,
  Title = '';
end

Blank = '';

for i = 5:length(rowlabels{1,1}),
  Blank = [Blank ' '];
end

length(Blank);

ColLabel = [];
for i=1:length(collabels),
  if ~isempty(collabels{i}),
    ColLabel{i} = collabels{i};
  end
end
[y,k] = sort(ColLabel);
ColLabel = ColLabel(k);

RowLabel = [];
for i=1:length(rowlabels),
  if ~isempty(rowlabels{i}),
    RowLabel{i} = rowlabels{i};
  end
end
[y,j] = sort(RowLabel);
RowLabel = RowLabel(j);

table = table(j,k);

T = sum(sum(table));                    % total count

% ------------------------------------- Raw count table -----------------
fprintf('%s - Raw counts\n',Title);
fprintf('         %s',Blank);
fprintf(' %4s ',ColLabel{:});
fprintf('RowCount Row%%ofTotal');
fprintf('\n');

for j=1:length(table(:,1)),
  fprintf('%5s     ', RowLabel{j});
  for k=1:length(table(1,:)),
    fprintf(' %4d ', table(j,k));
  end
  fprintf('  %5d', sum(table(j,:)));
  fprintf('  %5.1f%%\n', 100*sum(table(j,:))/T);
end

fprintf('Total    %s',Blank);
for k=1:length(table(1,:)),
  fprintf('%5.0f ', sum(table(:,k)));
end
fprintf('\n');

fprintf('Percent  ');
fprintf('%s',Blank);
for k=1:length(table(1,:)),
  fprintf('%5.0f%%', 100*sum(table(:,k))/T);
end
fprintf('\n%4d total count\n\n',T);

