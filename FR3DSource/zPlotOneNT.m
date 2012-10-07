function [void] = zPlotOneNT(NT,ViewParam)

X  = NT.Fit;

if isfield(ViewParam,'GlycoAtomSize'),
  GlycoAtomSize = ViewParam.GlycoAtomSize;
else
  GlycoAtomSize = 28;
end

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
  if strcmp(class(NT),'char'),
    LT = str2num(LT);
  end
else
  LT = 2.0;                              % doesn't work for some reason!
end

if isfield(ViewParam,'LabelBases'),
  LB = ViewParam.LabelBases;
else
  LB = 0;
end

if isfield(ViewParam,'Label'), 
  Label = ViewParam.Label;
elseif LB > 0,
  Label = [NT.Base NT.Number];
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

switch NT.Code
  case 1,
    col = [1 0 0];   % A is red 
  case 2,
    col = [1 0.8 0]; % C is yellow
  case 3,
    col = [0 1 0];   % G is green
  case 4, 
    col = [0 0 1];   % U is blue
  case 5,
    col = [155 48 255]/255; % modified are purple
  otherwise,
    col = [0 0 0];
end

gray = 0.5*[1 1 1];

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

switch NT.Code
case 1,
  i = [1 7 10 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [5 4 9 3 2 8];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);

  i = [7 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [15 6 14];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [9 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
case 2,
  i = [1 7 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [5 4 2 3];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);

  i = [7 10];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [8 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [12 6 13];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
case 3,
  i = [1 7 10 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [5 4 9 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [9 3 2 8];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);

  i = [7 13];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [4 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [16 11 15];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
case 4,
  i = [1 7 8 5 6];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [5 4 2 3];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);
  i = [1 2];
  plot3(X(i,1),X(i,2),X(i,3),'Color',col,'LineWidth',LT,'LineStyle',LS);

  i = [7 12];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [8 9];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
  i = [4 11];
  plot3(X(i,1),X(i,2),X(i,3),'Color',gray,'LineWidth',LT,'LineStyle',LS);
case 5,
  % we need a better algorithm to plot modified nucleotides!
  % we can use NT.Unit to tell which one it is, make a big switch statement
  % or read the CONECT information on each one in a database
%  i = [1:length(X(:,1)) 1];               % full circle

  D = zDistance(X,X);               % distances
  [i,j] = find((triu(D) <= 1.7) .* (triu(D) > 0));
  for k = 1:length(i),
    plot3([X(i(k),1) X(j(k),1)],[X(i(k),2) X(j(k),2)],[X(i(k),3) X(j(k),3)],'Color',col,'LineWidth',LT,'LineStyle',LS);
  end

  k = setdiff(1:length(X(:,1)),[i; j]);        % omitted nucleotides
  for kk = 1:length(k),
    i = k(kk);                                % disconnected atom
    [y,j] = sort(D(i,:));                     % j(2) is nearest atom
    plot3([X(i,1) X(j(2),1)],[X(i,2) X(j(2),2)],[X(i,3) X(j(2),3)],'Color',col,'LineWidth',LT,'LineStyle',LS);
  end

%  plot3(X(:,1),X(:,2),X(:,3),'.','Color',col,'MarkerSize',18);

end

if GlycoAtomSize > 0,
  if NT.Code <= 4,
    scatter3(X(1,1),X(1,2),X(1,3),GlycoAtomSize,col,'filled');% glycosidic atom
  end
end

if Sugar == 1,
 if NT.Code <= 4,                 % standard base
   Z = [NT.Sugar; NT.Fit(1,:)];
 elseif NT.Code == 5,             % modified base
   j = 1;
   for jj = length(NT.AtomName):-1:1,
     if strcmp(NT.AtomName{jj},'N1') || strcmp(NT.AtomName{jj},'N9'),
       j = jj;
     end
   end
   
%   D = zDistance(NT.Sugar,NT.Fit);
%   [i,j] = find(D == min(min(D)));% indices of nearest atoms

   Z = [NT.Sugar; NT.Fit(j,:)];
   if GlycoAtomSize > 0,
    scatter3(X(j,1),X(j,2),X(j,3),GlycoAtomSize,col,'filled');% glycosidic atom
   end
 else
   Z = [NT.Sugar; NT.Sugar(end,:)];
 end

 if length(NT.Sugar(:,1)) == 13,   % for some reason, some have 9

  k = [14 1 7 6 8 9 10 12]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  hold on
  k = [11 10]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [6 4 5]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [4 2 3]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [2 1]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);

  % connect backbone, if distance is small enough
  if norm(Z(13,1)-Z(10,1)) < 6,       % O3' from previous nucleotide is known
    k = [10 13];
    plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  end  

  if LSugar > 0,
    A = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P'};
    for j=1:12,
      text(Z(j,1)+0.1,Z(j,2),Z(j,3), A{j},'fontweight','bold','FontSize',LSugar);
    end
  end
 elseif length(NT.Sugar(:,1)) == 12,   % for some reason, some have 9

  k = [13 1 7 6 8 9 10 12]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  hold on
  k = [11 10]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [6 4 5]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [4 2 3]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);
  k = [2 1]; 
  plot3(Z(k,1),Z(k,2),Z(k,3),'Color',bc,'LineWidth',LT,'LineStyle',LS);

  if LSugar > 0,
    A = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P'};
    for j=1:12,
      text(Z(j,1)+0.1,Z(j,2),Z(j,3), A{j},'fontweight','bold','FontSize',LSugar);
    end
  end
 end
end

if LB > 0,
  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,Label,'fontweight','bold','FontSize',LB);
%  text(X(1,1)+0.5,X(1,2),X(1,3)+0.5,Label,'fontweight','bold','FontSize',LB);
end

if ShowBeta > 0 && isfield(NT,'Beta'),
  for j = 1:min(12,length(NT.Beta(:,1))),
    if NT.Beta(j,1) < Inf,
      text(NT.Sugar(j,1),NT.Sugar(j,2),NT.Sugar(j,3),num2str(round(NT.Beta(j,1))));
    end
  end
  for j = 13:length(NT.Beta),
    k = j - 12;
    if NT.Beta(k,1) < Inf,
      text(NT.Fit(j-12,1),NT.Fit(j-12,2),NT.Fit(j-12,3),num2str(round(NT.Beta(j,1))));
    end
  end
end
      

