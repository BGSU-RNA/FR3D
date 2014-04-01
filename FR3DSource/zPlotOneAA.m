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
  col = [255 65 150]/255;
end  




if isfield(ViewParam,'Color'),
  if length(ViewParam.Color) == 3,
    col = ViewParam.Color;
  end
end

if strcmp(LS,'-.'),
  col = 0.7*col;
  gray = 0.7*gray;
  LS = '-';
end

bc = gray;

hold on 

X  = AA.Loc;

plot3(X(:,1),X(:,2),X(:,3),'Color',col,'LineWidth',2,'LineStyle',LS);

if LB > 0,
  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
%  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
end

if ShowBeta > 0 && isfield(AA,'Beta'),
  for j = 1:length(AA.Beta(:,1)),
      text(AA.Loc(j,1),AA.Loc(j,2),AA.Loc(j,3),num2str(round(AA.Beta(j,1))));
  end
end
