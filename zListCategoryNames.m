% zListCategoryNames displays the category name for category Cat

Categories = [1 2 3 -3 4 -4 5 -5 6 -6 7 8 9 -9 10 -10 11 -11 12 -12 14  ...
21 -21 22 23 25 26 30 40:45 50 51];

for k=1:length(Categories),
  fprintf('%3d  %s\n',Categories(k), zCategoryName(Categories(k)));
end
