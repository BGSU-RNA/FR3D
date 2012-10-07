% zExamplarsAndInstances(File,Paircode,Cateogry) displays the best known
% representatives for interactions involving pairs with the given Paircode
% and interaction Category, together will all instances of that category

% Paircode and Category can be
% vectors.  It will loop through all possible Paircode, Category pairs from the
% two vectors.  If a certain category has subcategories, like 1, 1.1, and
% 1.2, it will loop through all of those.

% Here are some ways to run the program:

% zExemplarsAndInstances(File,7,21:23) paircode 7, cat's 21 to 23 (stacking)
% zExemplarsAndInstances(File,1:16,1)   all paircodes, category 1
% zExemplarsAndInstances(File,1:16,-12:23) all paircodes, all categories

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

function [void] = zExemplarsAndInstances(File,Paircode,Category)

% load exemplars -------------------------------------

  load('PairExemplars','Exemplar');

% specify parameters for viewing -------------------------------------------

  ViewParam.Mode      = 1; 
  ViewParam.Normal    = 1;
  ViewParam.ColorAxis = [-12 30];
  ViewParam.SortKeys  = [];
  ViewParam.Nearby    = 0;
  ViewParam.Sugar     = 1;
  ViewParam.ConnectSugar = 0;
  ViewParam.AtOrigin  = 1;
  ViewParam.Hydrogen  = 1;
  ViewParam.Sort      = 0;
  ViewParam.LineStyle = '-';

% loop through computer classifications ----------------------

for pc = 1:length(Paircode),
 for ca = 1:length(Category),
   for row = 1:length(Exemplar(:,Paircode(pc))),

     E = Exemplar(row,Paircode(pc));

     if ~isempty(E),
      if fix(E.Class) == Category(ca),
       fprintf('%s%s %s %s%s %s Category %3.1f\n',E.NT1.Base,E.NT1.Number,E.Pair.EdgeText,E.NT2.Base,E.NT2.Number,E.Filename,E.Class);
       zListPairData(E.Pair,1);

% display the exemplar pair ---------------------------------------------

       figure(1)
       cla
       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
       zPlotHydrogenBonds(E.NT1,E.NT2,E.Class,E.NT1.Rot,E.NT1.Fit(1,:));

       view(2)
       grid on
       axis equal

       switch fix(E.Class),
         case {1, 2},  axis([-2 10 -2 12 -5 5]);
         case 15,      axis([-5 3 -3 5 -5 5]);
         case 16,      axis([-4 4 -3 5 -3 5]);
         otherwise,    axis([-6 10 -6 10 -6 10]);
       end

       Title = [E.NT1.Base E.NT2.Base ' ' num2str(E.Class) ' ' strrep(E.Filename,'_','\_') ' '];
       Title = [Title E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
       title(Title);

       rotate3d on

% list pairs with the same interaction --------------------------

       P.Paircode  = Paircode(pc); 
       P.Category  = E.Class;
       P.Decimal   = 1;        % 1 - use 1.0 only; 0 - round to 1
       P.Group     = 1;        % hand, computer, both, etc.
       P.Sequential= 0;
       P.Expert    = 0;
       P.Context   = 3;

       VP.Mode      = 2;             
       VP.Color     = 12;
       VP.Normal    = 0;
       aa = Category(ca);
       ab = Category(ca) + sign(Category(ca));
       VP.ColorAxis = [min(aa,ab) max(aa,ab)];
       VP.Nearby    = 0;
       VP.Sugar     = 0;
       VP.Hydrogen  = 1;
       VP.az        = 51;
       VP.el        = 14;
       VP.LineStyle = '-';              % default - thick solid lines
       VP.Exemplars = 0;
       VP.ClassLimits = 2;
       VP.FigNum    = 2;

       VP.Sort      = 1;
       VP.SortKeys  = [20];
       VP.ListItems = [1 22 2 3 4 5 6 9 10 11 12 14 16 17];

       FR3D_PairViewer(File,P,VP);

% scatterplot of pairs with the same interaction --------------------------

       P.Paircode  = Paircode(pc); 
       P.Category  = fix(E.Class);
       P.Decimal   = 0;        % 1 - use 1.0 only; 0 - round to 1
       P.Group     = 1;        % hand, computer, both, etc.
       P.Sequential= 0;
       P.Expert    = 0;
       P.Inbox     = [-5 5 -5 5 -5 5 -1.1 1.1 -95 275];  % almost everything
       P.Context   = 3;

       VP.Mode      = 1;             
       VP.Color     = 12;
       VP.Normal    = 0;
       aa = Category(ca);
       ab = Category(ca) + sign(Category(ca));
       VP.ColorAxis = [min(aa,ab) max(aa,ab)];
       VP.Nearby    = 0;
       VP.Sugar     = 0;
       VP.Hydrogen  = 1;
       VP.az        = 51;
       VP.el        = 14;
       VP.LineStyle = '-';              % default - thick solid lines
       VP.Exemplars = 0;
       VP.ClassLimits = 2;
       VP.FigNum    = 2;

       VP.Sort      = 1;
       VP.SortKeys  = [20];
       VP.ListItems = [1 22 2 3 4 5 6 9 10 11 12 14 16 17];

       FR3D_PairViewer(File,P,VP);

       fprintf('This is the total number of matches over all subcategories.\n');

       figure(2)
       e = E.Pair.Displ;
       scatter3(e(1),e(2),e(3),28,'r','filled')
       title('Pairs that fall into this category, colored by subcategory, if any')
       figure(3)
       cut = 90;
       scatter3(mod(E.Pair.Ang+cut,360)-cut,E.Pair.Normal(3),E.Pair.Gap,28,'r','filled')
       hold on
       title('Pairs that fall into this category, colored by subcategory, if any')

% scatterplot of pairs which do not have the same interaction ------------

       P.Paircode  = Paircode(pc); 
       P.Category  = fix(E.Class);
       P.Decimal   = 0;        % 1 - use 1.0 only; 0 - round to 1
       P.Group     = 8;        % nearest exemplar matches
       P.Sequential= 0;
       P.Expert    = 0;
       P.Context   = 3;

       VP.Mode      = 1;             
       VP.Color     = 9;
       VP.Normal    = 0;
       aa = Category(ca);
       ab = Category(ca) + sign(Category(ca));
       VP.ColorAxis = [min(aa,ab) max(aa,ab)];
       VP.SortKeys  = [];
       VP.Nearby    = 0;
       VP.Sugar     = 0;
       VP.Hydrogen  = 1;
       VP.Sort      = 0;
       VP.az        = 51;
       VP.el        = 14;
       VP.LineStyle = '-';              % default - thick solid lines
       VP.Exemplars = 0;
       VP.ClassLimits = 2;
       VP.FigNum    = 4;

       FR3D_PairViewer(File,P,VP);

       fprintf('This is the total number of near matches over all subcategories.\n');
       figure(4)
       title('Pairs that do not fall into this category, colored by distance to exemplar')
       figure(5)
       title('Pairs that do not fall into this category, colored by distance to exemplar')

% view individual pairs that do not have the same classification -----------

disp('Press a key to start examining near misses for this (sub)category, closest to exemplar first')
pause

       P.Paircode  = Paircode(pc); 
       P.Category  = E.Class;
       P.Decimal   = 1;        % 1 - use 1.0 only; 0 - round to 1
       P.Group     = 8;        % nearest exemplar matches
       P.Sequential= 0;
       P.Expert    = 0;
       P.Context   = 3;

       VP.Mode      = 4;             
       VP.Color     = 2;
       VP.Normal    = 0;
       aa = Category(ca);
       ab = Category(ca) + sign(Category(ca));
       VP.ColorAxis = [min(aa,ab) max(aa,ab)];
       VP.Nearby    = 0;
       VP.Sugar     = 0;
       VP.Hydrogen  = 1;
       VP.az        = 51;
       VP.el        = 14;
       VP.LineStyle = '-';              % default - thick solid lines
       VP.Exemplars = 1;
       VP.ClassLimits = 1;
       VP.FigNum    = 6;

       VP.Sort      = 1;
       VP.SortKeys  = [20];
       VP.ListItems = [1 22 2 3 4 5 6 9 10 11 12 14 16 17];

       FR3D_PairViewer(File,P,VP);

% view individual pairs that have the same classification ------------------
disp('Press a key to start examining matches to this (sub)category, furthest from exemplar first')
pause

       P.Paircode  = Paircode(pc); 
       P.Category  = E.Class;
       P.Decimal   = 1;        % 1 - use 1.0 only; 0 - round to 1
       P.Group     = 1;        % nearest exemplar matches
       P.Sequential= 0;
       P.Expert    = 0;
       P.Context   = 3;

       VP.Mode      = 4;             
       VP.Color     = 2;
       VP.Normal    = 0;
       aa = Category(ca);
       ab = Category(ca) + sign(Category(ca));
       VP.ColorAxis = [min(aa,ab) max(aa,ab)];
       VP.Nearby    = 0;
       VP.Sugar     = 0;
       VP.Hydrogen  = 1;
       VP.az        = 51;
       VP.el        = 14;
       VP.LineStyle = '-';              % default - thick solid lines
       VP.Exemplars = 1;
       VP.ClassLimits = 1;
       VP.FigNum    = 6;

       VP.Sort      = 1;
       VP.SortKeys  = [-20];
       VP.ListItems = [1 22 2 3 4 5 6 9 10 11 12 14 16 17];

       FR3D_PairViewer(File,P,VP);


       fprintf('Press a key to go on to the next exemplar\n');
       pause

      end
     end
    end
  end
end
