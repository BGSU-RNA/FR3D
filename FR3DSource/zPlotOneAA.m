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

if isfield(ViewParam,'LabelAtoms'),
  LabelAtoms = ViewParam.LabelAtoms;
else
  LabelAtoms = 0;
end

if isfield(ViewParam,'ShowBeta'),
  ShowBeta = ViewParam.ShowBeta;
else
  ShowBeta = 0;
end

if isfield(ViewParam,'AADisplay'),
  AADisplay = ViewParam.AADisplay;
else
  AADisplay = 0;                          % 0 means show both groups, 1 means only One, 2 means only Two
end

gray = 0.5*[1 1 1];

col = [255 65 150]/255;
col = [	205  	41  	144]/255;             % default

if AADisplay == 0,
  GroupOneColor = col;
  GroupTwoColor = col;
  GroupThreeColor = col;
else
  GroupOneColor = [0 0 1];
  GroupTwoColor = gray;
  GroupThreeColor = col;
end

if AADisplay == 3,
  AADisplay = 0;                           % treat as 0 now that coloring is done
end

if isfield(ViewParam,'GroupOneColor') && length(ViewParam.GroupOneColor) == 3,
  GroupOneColor = ViewParam.GroupOneColor;
end

if isfield(ViewParam,'GroupTwoColor') && length(ViewParam.GroupTwoColor) == 3,
  GroupTwoColor = ViewParam.GroupTwoColor;
end

if isfield(ViewParam,'GroupThreeColor') && length(ViewParam.GroupThreeColor) == 3,
  GroupThreeColor = ViewParam.GroupThreeColor;
end

if strcmp(LS,'-.'),
%  col = 0.7*col;
  gray = 0.7*gray;
  LS = '-';
end

bc = gray;

hold on 

X  = AA.Loc;

if isfield(AA,'Unit') && isfield(AA,'Atom'),
  [ColorOne,ColorTwo,ColorThree] = zAAGroups(AA.Unit,AA.Atom);

  if any(AADisplay == [0 1]),
    k = ColorOne;
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupOneColor,'LineWidth',LT,'LineStyle',LS);
    if length(ColorOne) == 1 && length(ColorTwo) > 0,
      k = [ColorOne(end) ColorTwo(1)];
      plot3(X(k,1),X(k,2),X(k,3),'Color',GroupOneColor,'LineWidth',LT,'LineStyle',LS);
    end
  end

  if AADisplay == 0 && length(ColorOne) > 1 && length(ColorTwo) > 0,
    k = [ColorOne(end) ColorTwo(1)];
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupTwoColor,'LineWidth',LT,'LineStyle',LS);
  end

  if any(AADisplay == [0]) && length(ColorTwo) > 1,
    k = ColorTwo;
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupTwoColor,'LineWidth',LT,'LineStyle',LS);
  end

  if any(AADisplay == [0 2]) && length(ColorTwo) == 0 && length(ColorThree) == 1,
    k = [ColorOne(end) ColorThree(1)];
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupThreeColor,'LineWidth',LT,'LineStyle',LS);
  end

  if AADisplay == 0 && length(ColorTwo) > 0 && length(ColorThree) > 1,
    k = [ColorTwo(end) ColorThree(1)];
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupTwoColor,'LineWidth',LT,'LineStyle',LS);
  end

  if any(AADisplay == [0 2]) && length(ColorTwo) > 0 && length(ColorThree) == 1,
    k = [ColorTwo(end) ColorThree(1)];
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupThreeColor,'LineWidth',LT,'LineStyle',LS);
  end

  if any(AADisplay == [0 2]) && length(ColorThree) > 1,
    k = ColorThree;
    plot3(X(k,1),X(k,2),X(k,3),'Color',GroupThreeColor,'LineWidth',LT,'LineStyle',LS);
  end

  if LB > 0,
    text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
  %  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[AA.Unit AA.Number],'fontweight','bold','FontSize',LB);
  end

  if LabelAtoms > 0,
    for j = 1:length(X(:,1)),
      text(X(j,1),X(j,2),X(j,3),AA.Atom{j});
    end
  end

  if ShowBeta > 0 && isfield(AA,'Beta'),
    for j = 1:length(AA.Beta(:,1)),
        text(AA.Loc(j,1),AA.Loc(j,2),AA.Loc(j,3),num2str(round(AA.Beta(j,1))));
    end
  end
else
  D = zDistance(X,X);               % distances
  [i,j] = find((triu(D) <= 1.7) .* (triu(D) > 0));
  for k = 1:length(i),
    plot3([X(i(k),1) X(j(k),1)],[X(i(k),2) X(j(k),2)],[X(i(k),3) X(j(k),3)],'Color',col,'LineWidth',LT,'LineStyle',LS);
  end
end

