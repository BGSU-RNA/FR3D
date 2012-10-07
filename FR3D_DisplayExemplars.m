% zDisplayExamplars(Paircode,Cateogry) displays the best known
% representatives for interactions involving pairs with the given Paircode
% and interaction Category

% zDisplayExemplars(Paircode,Category), where Paircode and Category can be
% vectors, will loop through all possible Paircode, Category pairs from the
% two vectors.  If a certain category has subcategories, like 1, 1.1, and
% 1.2, it will loop through all of those.

% Here are some ways to run the program:

% zDisplayExemplars(7,15:18) paircode 7, categories 15 to 18 (stacking)
% zDisplayExemplars(1:16,1)   all paircodes, category 1
% zDisplayExemplars(1:16,-12:18) all paircodes, all categories

function [void] = zDisplayExemplars(Paircode,Category)

if nargin < 2
  Category = [-12:-1 1:23];
end

if nargin < 1
  Paircode = [1:16];
end

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

% loop through pairs and classifications ----------------------

for pc = 1:length(Paircode),
 for ca = 1:length(Category),
   for row = 1:length(Exemplar(:,Paircode(pc))),

     E = Exemplar(row,Paircode(pc));

     if fix(E.Class) == Category(ca),

       [Pair,s] = zClassifyPair(E.NT1,E.NT2);

fprintf('%10s %5s %s %5s   ',E.Filename, [E.NT1.Base E.NT1.Number],Pair.EdgeText,[E.NT2.Base E.NT2.Number]);

       figure(1)
       clf

       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
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

       fprintf('Press a key to go on\n');
       pause

      end
    end
  end
end
