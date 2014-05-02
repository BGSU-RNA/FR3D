function [void] = zPlotOneAA(AA,ViewParam)

if isfield(ViewParam,'LineStyle'),
  LS = ViewParam.LineStyle;
else
  LS = '-';
end

if isfield(ViewParam,'LineThickness'),   
  LT = ViewParam.LineThickness;
else
  LT = 2.0;                              % doesn't work for some reason!
end

if isfield(ViewParam,'LabelBases'),
  LB = ViewParam.LabelBases;
else
  LB = 0;
end

if isfield(ViewParam,'LabelSugar'),
  LSugar = ViewParam.LabelSugar;
else
  LSugar = 0;
end

if isfield(ViewParam,'ShowBeta'),
  ShowBeta = ViewParam.ShowBeta;
else
  ShowBeta = 0;
end

gray = 0.5*[1 1 1];

col = [255 65 150]/255;
col = [	205  	41  	144]/255;             % default

if strcmp(AA.Unit,'LYS'),
%  col = [255 65 150]/255;
end  

if isfield(ViewParam,'Color'),
  if length(ViewParam.Color) == 3,
%    col = ViewParam.Color;
  end
end

if strcmp(LS,'-.'),
%  col = 0.7*col;
  gray = 0.7*gray;
  LS = '-';
end

bc = gray;

hold on 

X  = AA.Loc;

D = zDistance(X,X);               % distances
[i,j] = find((triu(D) <= 1.7) .* (triu(D) > 0));
for k = 1:length(i),
  plot3([X(i(k),1) X(j(k),1)],[X(i(k),2) X(j(k),2)],[X(i(k),3) X(j(k),3)],'Color',col,'LineWidth',LT,'LineStyle',LS);
end

k = setdiff(1:length(X(:,1)),[i; j]);        % omitted atoms
for kk = 1:length(k),
  i = k(kk);                                % disconnected atom
  [y,j] = sort(D(i,:));                     % j(2) is nearest atom
  plot3([X(i,1) X(j(2),1)],[X(i,2) X(j(2),2)],[X(i,3) X(j(2),3)],'Color',col,'LineWidth',LT,'LineStyle',LS);
end

if LB > 0,
  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
%  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
end

if ShowBeta > 0 && isfield(AA,'Beta'),
  for j = 1:length(AA.Beta(:,1)),
      text(AA.Loc(j,1),AA.Loc(j,2),AA.Loc(j,3),num2str(round(AA.Beta(j,1))));
  end
end
