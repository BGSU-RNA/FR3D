% zExemplarTable(Cateogry) displays the best known representatives for interactions involving all pairs in interaction category(ies) Category
% Example:  zExemplarTable(1,3.5)
% Example:  zExemplarTable(1:12,3.5)

% Coarse = 1 produces a 4-color version of the isodiscrepancy figure
% Subcat = 1, include subcategories, Subcat = 0, don't

function [T, U] = zExemplarTable(Category,Coarse,Subcat,Verbose)

if nargin < 1,
  Category = 1;
end

if nargin < 2,
  Coarse = 1;
end

if nargin < 3,
  Subcat = 1;
end

if nargin < 4,
  Verbose = 2;
end

% load exemplars -------------------------------------

load('PairExemplars','Exemplar');

% loop through computer classifications, accumulating exemplars ------------

B(1) = Exemplar(1,1);
B(1).HydrogenClass = 0;
B(1).subplot  = 0;
B(1).original = 0;

Lab = [];
Cat = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

for c = 1:length(Category),
for c1 = 1:4,
 for c2 = 1:4,
  pc  = 4*(c2-1)+c1;                       % current paircode
  for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar

    E = Exemplar(r,pc);                    % current exemplar

    if ~isempty(E.NT1),                    % non-empty entry of Exemplar
    if  any(abs(E.Class) == Category(c)) || ...
       (any(fix(abs(E.Class)) == Category(c)) && (Subcat == 1)),
    if (E.Count >= 0),
    
      [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);

      % ------- in some symmetric families, produce the symmetric version
      if (any(pc == [5 7 9 13 14 15])) && any(fix(E.Class) == [1 2 7 8]),
        E.Class = -E.Class;
        [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);
        E.Class = -E.Class;
      end

      % ------- in some symmetric families, store AA, CC, GG, UU pairs twice

      if (E.NT1.Code == E.NT2.Code) && any(fix(E.Class) == [1 7 8]),
        E.HydrogenClass = 0;
        E.subplot  = 0;
        E.original = 0;
        E.Class = -E.Class;
        [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);
      end


    end
    end
    end
  end
 end
end
end

% remove first exemplar --------------------

B   = B(2:end);
Lab = Lab(2:end);
Cat = Cat(2:end);

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
ViewParam.LabelBases= 8;                    % font size

% -------------------------------------- Plot exemplars

if Verbose > 1,
  for ca = 1:length(Category),
    figure(fix(Category(ca)))
    clf
    plotted = zeros(16,1);                   % keep track of which are plotted
    for m = 1:length(B),
     if (B(m).original == 1) && ((abs(B(m).Class) == Category(ca)) ...
        || ((abs(fix(B(m).Class)) == Category(ca)) && (Subcat == 1))),

        pc2 = B(m).subplot;
        E   = B(m);
        % display the exemplar pair -----------------------------------------

       if abs(E.Class - fix(E.Class)) == 0,
         ViewParam.LineStyle = '-';
       elseif abs(E.Class - fix(E.Class)) > 0.29,
         ViewParam.LineStyle = '.';
       elseif abs(E.Class - fix(E.Class)) > 0.19,
         ViewParam.LineStyle = ':';
       elseif abs(E.Class - fix(E.Class)) > 0.09,
         ViewParam.LineStyle = '--';
       end

       subplot(4,4,pc2);

       if plotted(pc2) > 0,
         ViewParam.LabelBases = 0;                  % don't display now
         xlab{pc2} = [xlab{pc2} ', ' num2str(E.Count)];  % append subcat count
       else,
         ViewParam.LabelBases = 8;                  % font size
         xlab{pc2} = ['Count: ' num2str(E.Count)];
         Title{pc2} = [E.NT1.Base E.NT2.Base zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ' ' strrep(E.Filename,'_','\_') ' '];
         Title{pc2} = [Title{pc2} E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
         CP = norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:));     % c1'-c1' dist
         Title{pc2} = [Title{pc2} ' ' num2str(CP)];
       end

       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
       zPlotHydrogenBonds(E.NT1,E.NT2,E.HydrogenClass,E.NT1.Rot,E.NT1.Fit(1,:));

       view(2)
       grid off
       axis equal
       %axis tight
       axis fill
       xlabel(xlab{pc2});
       title(Title{pc2});

%       a = axis;
%       text(0.3*a(1)+0.7*a(2),0.8*a(3)+0.2*a(4),['Count: ' num2str(E.Count)]);

       rotate3d on
       plotted(pc2) = 1;

     end
    end
  end

 h = gcf;
 orient landscape
% print(h,'-depsc2',['Isostericity\Exemplars' zEdgeText(Category(ca))]);

end

% -------------------------------------- Compare basepairs against each other

if exist('B'),
D = [];
RD = [];
G = zeros((length(B)^2-length(B))/2,7);
j = 1;

for a = 1:length(B),
  for b = a:length(B),

    % calculate isodiscrepancy

    D(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT1,B(b).NT2);

if B(a).Class == B(b).Class,
%  LookupTable{B(a).Class}(B(a).NT1.Code,B(a).NT2.Code,B(b).NT1.Code,B(b).NT2.Code) = D(a,b);
%  LookupTable{B(a).Class}(B(a).NT2.Code,B(a).NT1.Code,B(b).NT2.Code,B(b).NT1.Code) = D(a,b);
end
  

    D(b,a) = D(a,b);

    RD(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT2,B(b).NT1);
    RD(b,a) = RD(a,b);               % reverse order, for family distance

    j = j + 1;

  end
end

% ------------------------------------------ Print table of isodiscrepancies
% ------------------------------------------ Display graph of isodiscrepancies

if length(Category) > 1,
  p = zClusterGraph(D, Lab, 12, 1:length(B), 0);
else
  p = zClusterGraph(D, Lab, [12 2], [], Verbose);
  if length(Category) == 1 && Category(1) == 1,
    p = [9 18  4  16   8   7  14  13  15   3  17  10   12  11 1  6   5   2];
  end
  zClusterGraph(D, Lab, [12 2], p, Verbose);
end



colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 8]);
colorbar('location','eastoutside');

Title = [];
for i = 1:length(Category),
  Title = [Title ' ' zEdgeText(abs(Category(i)))];
end
Title = [Title 'family isodiscrepancy map'];
if Subcat > 0,
  Title = [Title ' with subcategories'];
end
title(Title,'FontSize',10);

saveas(gcf,['Isostericity' filesep Title '.pdf'],'pdf');
saveas(gcf,['Isostericity' filesep Title '.png'],'png');

if Coarse > 0,
  figure
  w = [1.55 2.25 3.3 4.5 5.5 8];
  w = [2 2.6 3.3 4.5 5.5 8];
  E = (0 * (D < w(1))) + w(2) * (D > w(1)) .* (D < w(3)) + w(4) * (D > w(3)) .* (D < w(5)) + w(6) * (D > w(5));
  zClusterGraph(E,Lab,[12 2],p,0);

  colormap('default');
  map = colormap;
  map = map((end-8):-1:8,:);
  colormap(map);
  caxis([0 8]);
  colorbar('location','eastoutside');

  Title = [];
  for i = 1:length(Category),
    Title = [Title ' ' zEdgeText(abs(Category(i)))];
  end
  Title = [Title 'family isodiscrepancy map'];
  if Subcat > 0,
    Title = [Title ' with subcategories'];
  end
  title(Title,'FontSize',10);
end

% ------------------------------------------ Provide cell output of table

%AAcHS  I1/I2    9.0 1GRZ   32
%123456789012345678901234567890

T = {};
T{1,1} = zEdgeText(Category(1));
T{1,1} = T{1,1}(1:3);
T{1,2} = 'Family';
T{1,3} = 'LSW 2000 subfamily';
T{1,4} = 'Structure';
T{1,5} = 'Count';
for i = 1:length(Lab),
  T{1,5+i} = Lab{p(i)}(1:2);
  T{i+1,1} = Lab{p(i)}(1:2);
  T{i+1,2} = Lab{p(i)}(3:5);
  T{i+1,3} = Lab{p(i)}(8:12);
  T{i+1,4} = Lab{p(i)}(21:24);
  T{i+1,5} = str2num(Lab{p(i)}(26:29));  
  for j = 1:length(Lab),
    T{i+1,5+j} = D(p(i),p(j));
  end
end

% ------------------------------------------ Provide cell output of 4x4 tables

U = {};
[y,i] = sort(Lab);
SLab = Lab(i);
SLab'
DD = D(i,i);                          % re-order alphabetically
SLab{end+1} = '  ';

r = 1;
k = 1;
while k <= length(Lab),
  U{r,1} = [SLab{k}(1:2) ' ' SLab{k}(3) upper(SLab{k}(4:5))];
  U{r,2} = 'A';
  U{r,3} = 'C';
  U{r,4} = 'G';
  U{r,5} = 'U';
  U{r+1,1} = 'A';
  U{r+2,1} = 'C';
  U{r+3,1} = 'G';
  U{r+4,1} = 'U';

  if strcmp(SLab{k}(1:2),SLab{k+1}(1:2)),       % AA, CC, GG, or UU
    for j = 1:length(Lab),
      a = find(SLab{j}(1) == 'ACGU');
      b = find(SLab{j}(2) == 'ACGU');
      if isempty(U{r+a,1+b}),
        U{r+a,1+b} = min(DD(k,j),DD(k+1,j));     % min of this row and next
      else
        U{r+a,1+b} = min([U{r+a,1+b} DD(k,j) DD(k+1,j)]);
      end
    end
    k = k + 1;
  else
    for j = 1:length(Lab),
      a = find(SLab{j}(1) == 'ACGU');
      b = find(SLab{j}(2) == 'ACGU');
      if isempty(U{r+a,1+b}),
        U{r+a,1+b} = DD(k,j);
      else
        U{r+a,1+b} = min(U{r+a,1+b},DD(k,j));
      end
    end
  end
  k = k + 1;
  r = r + 6;                                     % current row
end

% ----------------------------------------- Calculate avg dist btw families

if length(Category) > 1,

Fam{1}  = 'cWW';
Fam{2}  = 'tWW';
Fam{3}  = 'cWH';
Fam{4}  = 'tWH';
Fam{5}  = 'cWS';
Fam{6}  = 'tWS';
Fam{7}  = 'cHH';
Fam{8}  = 'tHH';
Fam{9}  = 'cHS';
Fam{10} = 'tHS';
Fam{11} = 'cSS';
Fam{12} = 'tSS';

MeanD = [];

minD = min(D,RD);

for p = 1:12,
  i = find(ismember(Cat,upper(Fam{p})));
  if ~isempty(i),
    for q = 1:12,
      j = find(ismember(Cat,upper(Fam{q})));
      if ~isempty(j),
%        MeanD(p,q) = mean(mean(minD(i,j)));
        MeanD2(p,q) = mean(mean(D(i,j)));
        MeanD3(p,q) = mean(mean(RD(i,j)));
%        MeanD2(p,q) = median(median(D(i,j)));
%        MeanD3(p,q) = median(median(RD(i,j)));
        MeanD(p,q)  = min(MeanD2(p,q),MeanD3(p,q));
      end
    end
  end
end

h = find(sum(MeanD) > 0);                 % which pairs are actually present

figure

% Neocles says:
% Separate the 2 groups: Group I: cWW, cWS, tWH, cHH, tHS, cSS (1 5 4 7 10 11)
% Group II: tWW, tWS, cWH, tHH, cHS, tSS, (2 6 3 8 9 12)

p = [7 4 1 10 5 11  3 9 6 12 2 8];
p = zClusterGraph(MeanD(h,h),Fam(h),[3 3],p,2);

%p = zClusterGraph(MeanD(h,h),Fam(h),[3 3]);
%p

colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 12]);
colorbar('location','eastoutside');
title('Average IsoDiscrepancies between geometric families');

end
end

% --------------------------------------------------------------------------

function [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat)

m = length(B) + 1;

pc = 4*(E.NT2.Code-1)+E.NT1.Code;                  % paircode

E.HydrogenClass = E.Class;

if (E.Class < 0),
  T     = E.NT1;
  E.NT1 = E.NT2;
  E.NT2 = T;
  if E.NT1.Code ~= E.NT2.Code,
    E.Class = -E.Class;
  end
end

if (pc == 7),                                      % reverse GC pairs
  T     = E.NT1;
  E.NT1 = E.NT2;
  E.NT2 = T;
end

pc2 = 4*(E.NT1.Code-1)+E.NT2.Code;                 % paircode when reversed

E.subplot = pc2;
if ~isfield(E,'original'),
  E.original = 1;
end

ic = zIsostericSubgroups(E.NT1.Code,E.NT2.Code,abs(fix(E.Class)));

B(m) = E;                    % store this exemplar for isodisc calc
Lab{m} = [E.NT1.Base E.NT2.Base zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ' ' ic sprintf('  %5.1f',E.Class) ' ' E.Filename ' ' sprintf('%4d',E.Count)];
Cat{m} = upper(zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code));


return  

zExemplarTable(1,3.5,0,1);
zExemplarTable(2,3.5,0,1);
zExemplarTable(3,3.5,0,1);
zExemplarTable(4,3.5,0,1);
zExemplarTable(5,3.5,0,1);
zExemplarTable(6,3.5,0,1);
zExemplarTable(7,3.5,0,1);
zExemplarTable(8,3.5,0,1);
zExemplarTable(9,3.5,0,1);
zExemplarTable(10,3.5,0,1);
zExemplarTable(11,3.5,0,1);
zExemplarTable(12,3.5,0,1);
zExemplarTable(13,3.5,0,1);
zExemplarTable([1:12],3.5,0,0);
zExemplarTable([1 5 4 7 10 11 2 6 3 8 9 12],3.5,0,0);

zExemplarTable([1 5],3.5,0);
zExemplarTable(1:6,3.5,0);
zExemplarTable(7:12,3.5,0);

% Neocles says:
% Separate the 2 groups: Group I: cWW, cWS, tWH, cHH, tHS, cSS
% Group II: tWW, tWS, cWH, tHH, cHS, tSS, 

