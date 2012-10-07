function [void] = zPlotOneNT(NT,ViewParam)

X  = NT.Fit;

if isfield(ViewParam,'Sugar'),
  Sugar = ViewParam.Sugar;
else
  Sugar = 0;
end

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

switch NT.Code
  case 1,
    col = [1 0 0];   % A is red 
  case 2,
    col = [1 0.8 0]; % C is yellow
  case 3,
    col = [0 1 0];   % G is green
  case 4, 
    col = [0 0 1];   % U is blue
  otherwise,
    col = [0 0 0];
end

gray = 0.5*[1 1 1];

if isfield(ViewParam,'Color'),
  if length(ViewParam.Color) == 3,
    col = ViewParam.Color;
  end
end

hold on 

switch NT.Code
case 1,
  i = [13 1 7 10 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [5 4 9 3 2 8];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);

  i = [7 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [15 6 14];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [9 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
case 2,
  i = [9 1 7 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [5 4 2 3];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);

  i = [7 10];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [8 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [12 6 13];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
case 3,
  i = [14 1 7 10 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [5 4 9 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [9 3 2 8];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);

  i = [7 13];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [4 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [16 11 15];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
case 4,
  i = [10 1 7 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [5 4 2 3];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',2,'LineStyle',LS);

  i = [7 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [8 9];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
  i = [4 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',2,'LineStyle',LS);
end

scatter3(X(1,1),X(1,2),X(1,3),28,col,'filled');   % glycosidic atom

if Sugar == 1,

if length(NT.Sugar(:,1)) == 12,   % for some reason, some have 9
  Z = [NT.Sugar; NT.Fit(1,:)];

  k = [13 1 7 6 8 9 10 12]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'k','LineWidth',2,'LineStyle',LS);
  hold on
  k = [11 10]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'k','LineWidth',2,'LineStyle',LS);
  k = [6 4 5]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'k','LineWidth',2,'LineStyle',LS);
  k = [4 2 3]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'k','LineWidth',2,'LineStyle',LS);
  k = [2 1]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'k','LineWidth',2,'LineStyle',LS);
end

end

if LB > 0,
  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[NT.Base NT.Number],'fontweight','bold','FontSize',LB);
%  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,[NT.Base NT.Number],'fontweight','bold','FontSize',LB);
end
