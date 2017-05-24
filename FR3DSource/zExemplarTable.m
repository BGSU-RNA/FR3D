% zExemplarTable(Cateogry) displays the best known representatives for interactions involving all pairs in interaction category(ies) Category

% zExemplarTable(1,0,0,1)

% Coarse  = 1 produces a 4-color version of the isodiscrepancy figure
% Subcat  = 1, include subcategories, Subcat = 0, don't
% Verbose = 1, shows some things
% Verbose = 2, shows 4x4 plots of all instances in the category

% T is the full IDI map for the category
% U is a set of 4x4 tables

function [T, U, FNS, CSV] = zExemplarTable(Category,Coarse,Subcat,Verbose)

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

load([pwd filesep 'FR3DSource' filesep 'PairExemplars'],'Exemplar');

% loop through computer classifications, accumulating exemplars ------------

B(1) = Exemplar(1,1);
B(1).HydrogenClass = 0;
B(1).subplot  = 0;
B(1).original = 0;

Lab = [];
Cat = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

for c = 1:length(Category),                % requested categories
  for c1 = 1:4,                              % loop through code of base 1
   for c2 = 1:4,                             % loop through code of base 2
    pc  = 4*(c2-1)+c1;                       % current paircode
    for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar

      E = Exemplar(r,pc);                    % current exemplar

      if ~isfield(E,'Source'),
        fprintf('No source field for this exemplar:\n');
        E
      end

      if ~isempty(E.NT1),                    % non-empty entry of Exemplar
        Epc = 4*(E.NT2.Code-1)+E.NT1.Code;   % paircode of exemplar itself

        if Epc ~= pc,                        % should not happen
          fprintf('Paircode disagreement, %d versus %d in row %d\n', pc, Epc, r);
        end

        if  any(abs(E.Class) == Category(c)) || ...
           (any(fix(abs(E.Class)) == Category(c)) && (Subcat == 1)),
          if (E.Count >= 0) && Epc == pc,

            [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);

            % ------- in some symmetric families, produce the symmetric version
            if (any(pc == [5 7 9 13 14 15])) && any(fix(E.Class) == [1 2 7 8 14]),
              E.Class = -E.Class;
              [B,Lab,Cat] = AddExemplar(E,B,Lab,Cat,Subcat);
              E.Class = -E.Class;
            end

            % ------- in some symmetric families, store AA, CC, GG, UU pairs twice

            if (E.NT1.Code == E.NT2.Code) && any(fix(E.Class) == [1 7 8 14]),
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

% remove first exemplar, just a placeholder --------------------

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
    figure(30+fix(Category(ca)))             % clear figures now, alternate below
    clf
    figure(fix(Category(ca)))                % for plotting basepars
    clf
    plotted = zeros(16,1);                   % keep track of which are plotted
    for m = 1:length(B),                     % go through all stored exemplars
      if (B(m).original == 1) && ((abs(B(m).Class) == Category(ca)) ...
        || ((abs(fix(B(m).Class)) == Category(ca)) && (Subcat == 1))),

        figure(fix(Category(ca)))

        pc2 = B(m).subplot;                  % which subplot this goes in
        E   = B(m);
        % display the exemplar pair -----------------------------------------

        ViewParam.Sugar = 0;

        if abs(E.Class - fix(E.Class)) == 0,
          ViewParam.LineStyle = '-';
          ViewParam.Sugar = 1;
        elseif abs(E.Class - fix(E.Class)) > 0.29,
          ViewParam.LineStyle = '.';
        elseif abs(E.Class - fix(E.Class)) > 0.19,
          ViewParam.LineStyle = ':';
        elseif abs(E.Class - fix(E.Class)) > 0.09,
          ViewParam.LineStyle = '--';
        end

        ViewParam.LineStyle = '-';             % make them all the same

        ha = subplot(4,4,pc2);            % tell where to plot this
        p = get(ha, 'pos');
        cp = 0.01;                             % change in percentage; don't go too big
        p(1) = p(1) - cp;
        p(2) = p(2) - cp;
        p(3) = p(3) + 2*cp;
        p(4) = p(4) + 2*cp;
        set(ha, 'pos', p);

        axis off

        if plotted(pc2) > 0,
          ViewParam.LabelBases = 0;                       % don't display now
          xlab = [xlab ', ' num2str(E.Count)];  % append subcat count
          Title{pc2} = [Title{pc2} ', ' num2str(E.Count)];
        else
          ViewParam.LabelBases = 8;                       % font size
          xlab = ['Count: ' num2str(E.Count)];

          g = strfind(E.NT1.ID,'|');
          if length(g) < 3,
            FN = E.Filename;
          else
            FN = E.NT1.ID(1:(g(3)-1));
          end
          FN = strrep(FN,'CURATED_','');
          FN = strrep(FN,'_EXEMPLAR','');
          FN = strrep(FN,'MODE','Model');

          Title{pc2} = [E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number ' ' zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ' ' FN ' '];
          Title{pc2} = [Title{pc2} ' ' sprintf('%d',E.Count)];
          Title{pc2} = strrep(Title{pc2},'  ',' ');
        end

        F.NT(1) = E.NT1;
        F.NT(2) = E.NT2;
        F.Filename = E.Filename;
        zDisplayNT(F,[1 2],ViewParam);
        hold on
        zPlotHydrogenBonds(E.NT1,E.NT2,E.HydrogenClass,E.NT1.Rot,E.NT1.Fit(1,:));

        view(2)
        grid off
        axis equal
        %axis tight
        axis fill

        ax = axis;
%        text((ax(1)+ax(2))/2,ax(3)-0.1*(ax(4)-ax(3)),xlab,'fontsize',8,'horizontalalignment','center')
%       xlabel(xlab);
%       ylabel(Title{pc2},'fontsize',7);
%        title('');
        title(Title{pc2},'fontsize',7);
%        rotate3d on
        plotted(pc2) = 1;

       % ------------------------------ plot and save each basepair separately

       figure(55)
       clf
       VP = ViewParam;
       VP.LineThickness = 4;
       VP.LabelBases = 0;
       zDisplayNT(F,[1 2],VP);
       hold on
       title('');
       zPlotHydrogenBonds(E.NT1,E.NT2,E.HydrogenClass,E.NT1.Rot,E.NT1.Fit(1,:));
       view(2)
       grid off
       axis off
       axis equal
       %axis tight
       axis fill
       FN = ['Exemplars' filesep zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) '_' E.NT1.Base E.NT2.Base '_exemplar.png'];
       saveas(gcf,FN,'png');

%       [X,map] = imread(FN);
%       Y = X(67:810,150:1100,:);
%       imwrite(Y,FN);

       % ------------------------------ plot glycosidic bonds

       figure(30+fix(Category(ca)))

       c1 = F.NT(1).Code;
       c2 = F.NT(2).Code;
       pc  = 4*(c2-1)+c1;                       % current paircode

       Shifts = zeros(16,2);

       if fix(Category(ca)) == 1,
         Shifts(7,2) = 0.2;
         Shifts(4,2) = -0.1;
         Shifts(10,2) = -0.4;
       end

       Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

       L1 = Lim(2,F.NT(1).Code);
       L2 = Lim(2,F.NT(2).Code);

       if (F.NT(1).Code == 4 && F.NT(2).Code == 4 && Category(ca) == 1) || ...
          (F.NT(1).Code == 4 && F.NT(2).Code == 1 && Category(ca) == 4),
         V = ViewParam;
         V.Sugar = 0;
         V.LabelBases = 0;
         V.GlycoAtomSize = 0;
         zDisplayNT(F,[1 2],V)
       end
       hold on

       R = F.NT(1).Rot;             % Rotation matrix for first base
       S = F.NT(1).Fit(1,:);        % Location of glycosidic atom
       NT1Fit = (F.NT(1).Fit   - ones(L1,1)*S) * R; % rotated into position
       NT1Sug = (F.NT(1).Sugar(1,:) - S) * R; % rotated into position
       NT2Fit = (F.NT(2).Fit   - ones(L2,1)*S) * R; % rotated into position
       NT2Sug = (F.NT(2).Sugar(1,:) - S) * R; % rotated into position

       H  = [13 9 14 10];                     % row of Fit containing H1/H9
       h  = H(F.NT(1).Code);
       h2 = H(F.NT(2).Code);

       plot([NT2Fit(1,1) NT2Fit(h2,1)],[NT2Fit(1,2) NT2Fit(h2,2)],'r', 'LineWidth', 2);

       if Category(ca) == 1 && ((c1 == 4 && c2 == 1) || (c1 == 2 && c2 == 3)),
         if c1 < c2,
           text(NT2Fit(h2,1)+Shifts(pc,1),NT2Fit(h2,2)+Shifts(pc,2),[' UA, CG']);
         end
       elseif Category(ca) == 1 && ((c1 == 3 && c2 == 1) || (c1 == 1 && c2 == 3)),
         if c1 < c2,
           text(NT2Fit(h2,1)+Shifts(pc,1),NT2Fit(h2,2)+Shifts(pc,2),[' GA, AG']);
         end
       else
         text(NT2Fit(h2,1)+Shifts(pc,1),NT2Fit(h2,2)+Shifts(pc,2),[' ' F.NT(1).Base F.NT(2).Base]);
       end

       plot([NT1Fit(1,1) NT1Fit(h,1)],[NT1Fit(1,2) NT1Fit(h,2)],'k', 'LineWidth', 2);
       plot(NT2Fit(h2,1), NT2Fit(h2,2), '.r', 'MarkerSize', 16);
       hold on
       plot(NT1Fit(h,1), NT1Fit(h,2), '.k', 'MarkerSize', 16);
       text(NT1Fit(h,1), NT1Fit(h,2), ' C1''');

       if F.NT(1).Code == F.NT(2).Code && any(Category(ca) == [1 2 7 8]),
         F.NT = F.NT([2 1]);          % change order of nucleotides
         R = F.NT(1).Rot;             % Rotation matrix for first base
         S = F.NT(1).Fit(1,:);        % Location of glycosidic atom
         NT1Fit = (F.NT(1).Fit   - ones(L1,1)*S) * R; % rotated into position
         NT1Sug = (F.NT(1).Sugar(1,:) - S) * R; % rotated into position
         NT2Fit = (F.NT(2).Fit   - ones(L2,1)*S) * R; % rotated into position
         NT2Sug = (F.NT(2).Sugar(1,:) - S) * R; % rotated into position

         H  = [13 9 14 10];
         h  = H(F.NT(1).Code);
         h2 = H(F.NT(2).Code);

         plot([NT2Fit(1,1) NT2Fit(h2,1)],[NT2Fit(1,2) NT2Fit(h2,2)],'r', 'LineWidth', 2);
         plot(NT2Fit(h2,1), NT2Fit(h2,2), '.r', 'MarkerSize', 16);
         text(NT2Fit(h2,1)+Shifts(pc,1),NT2Fit(h2,2)+Shifts(pc,2),[' ' F.NT(1).Base F.NT(2).Base]);
       end

       title('');
       axis equal

     end
    end

    % -------------------------------------- Save images of basepairs

    figure(fix(Category(ca)))
    orient landscape
    saveas(gcf,['Isostericity' filesep zEdgeText(Category) '_basepairs.png'],'png')
    saveas(gcf,['Isostericity' filesep zEdgeText(Category) '_basepairs.pdf'],'pdf')
    FNS{1} = ['Isostericity' filesep zEdgeText(Category) '_basepairs.pdf'];

    % -------------------------------------- Save images of glycosidic bonds

    figure(30+fix(Category(ca)))
    view(2)
    grid off
    axis off
    orient landscape
    title([zEdgeText(Category) 'glycosidic bonds, relative to glycosidic bond of first base']);
    saveas(gcf,['Isostericity' filesep zEdgeText(Category) 'PairGlydosidicBonds.png'],'png')
    saveas(gcf,['Isostericity' filesep zEdgeText(Category) 'PairGlydosidicBonds.pdf'],'pdf')
    FNS{2} = ['Isostericity' filesep zEdgeText(Category) 'PairGlydosidicBonds.pdf'];

  end
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
    D2(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT1,B(b).NT2,1);

    D(b,a) = D(a,b);
    D2(b,a) = D2(a,b);

    RD(a,b) = zIsoDiscrepancy(B(a).NT1,B(a).NT2,B(b).NT2,B(b).NT1);
    RD(b,a) = RD(a,b);               % reverse order, for family distance

    j = j + 1;

  end
end

% ------------------------------------------ Print table of isodiscrepancies
% ------------------------------------------ Display graph of isodiscrepancies

figure(20)

D = D2;                                    % use idealized glycosidic bond!

Lab
Lab'

if length(Category) > 1,
  p = zClusterGraph(D, Lab, 24, 1:length(B), 0);
else
  p = zClusterGraph(D, Lab, [24 2], [], Verbose);
  if length(Category) == 1 && Category(1) == 1 && length(D(:,1)) == 18,
    p = [9 18  4  16   8   7  14  13  15   3  17  10   12  11 1  6   5   2];
  end
  zClusterGraph(D, Lab, [24 2], p, Verbose);
end

colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 8]);
colorbar('location','eastoutside');

Title = [];
for i = 1:length(Category),
  Title = [Title zEdgeText(abs(Category(i)))];
end
Title = [Title 'family isodiscrepancy map'];
if Subcat > 0,
  Title = [Title ' with subcategories'];
end
title(Title,'FontSize',12);

% Title = Title(2:end);

orient landscape
saveas(gcf,['Isostericity' filesep Title '.pdf'],'pdf');
saveas(gcf,['Isostericity' filesep Title '.png'],'png');
FNS{3} = ['Isostericity' filesep Title '.pdf'];

if Coarse > 0,
  figure
  w = [1.55 2.25 3.3 4.5 5.5 8];
  w = [2 2.6 3.3 4.5 5.5 8];
  E = (0 * (D < w(1))) + w(2) * (D >= w(1)) .* (D < w(3)) + w(4) * (D >= w(3)) .* (D < w(5)) + w(6) * (D >= w(5));
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
  title(Title,'FontSize',12);
end

% ------------------------------------------ Provide cell output of table

%AAcHS  I1/I2    9.0 1GRZ   32
%123456789012345678901234567890

T = {};
T{1,1} = 'Family';
T{1,2} = 'Base combination';
T{1,3} = 'FragmentFilename';
T{1,4} = 'PDBID';
T{1,5} = 'NT1ID';
T{1,6} = 'NT2ID';
T{1,7} = 'Count';
T{1,8} = 'LW 2002 subgroup';
T{1,9} = 'Resolution';
T{1,10} = 'C1*-C1* distance';

% Family,Base combination,FragmentFilename,PDBID,NT1ID,NT2ID,Count,LW 2002 subgroup,Resolution,C1*-C1* distance,IDI,IDI,IDI,...

Catt  = [num2str(abs(Category(1))) '.'];

for i = 1:length(Lab),
  T{1,10+i} = [B(p(i)).NT1.Base B(p(i)).NT2.Base];
end

for i = 1:length(Lab),
  E = B(p(i));

  T{i+1,1} = strrep(zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code),' ','');
  T{i+1,2} = [E.NT1.Base E.NT2.Base];
  if E.NT1.Code == E.NT2.Code && E.Class < 0,
    T{i+1,3} = [T{i+1,1} '_' T{i+1,2} '_2_Exemplar.pdb'];
  else
    T{i+1,3} = [T{i+1,1} '_' T{i+1,2} '_Exemplar.pdb'];
  end
  T{i+1,4} = E.Filename;
  T{i+1,5} = E.NT1.ID;
  T{i+1,6} = E.NT2.ID;
  T{i+1,7} = num2str(E.Count);

  Cate = strrep(Lab{p(i)}(8:12),'I',Catt);
  Cate = strrep(Cate,'i',['i' Catt]);
  Cate = strrep(Cate,'(','[');
  Cate = strrep(Cate,')',']');
  T{i+1,8} = strrep(Cate,' ','');                 % isosteric subgroup
  T{i+1,9} = num2str(E.Resolution);
  T{i+1,10} = sprintf('%0.6f',norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:)));

  for j = 1:length(Lab),
    T{i+1,10+j} = sprintf('%0.6f',D(p(i),p(j)));
  end
end

[s,t] = size(T);
for i = 1:s,
  for j = 1:t,
    if j > 1,
      CSV{i} = [CSV{i} ',' T{i,j}];
    else
      CSV{i} = [T{i,j}];
    end
  end
end

% ------------------------------------------ Write out CSV file

fid = fopen(['Exemplars' filesep zEdgeText(Category(1)) '_data.txt'],'w');
for i = 1:length(CSV),
  fprintf(fid,'%s\n',CSV{i});
end
fclose(fid)

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

size(Cat)
size(h)
size(MeanD)

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

m = length(B) + 1;                            % add to the end of the list

E.HydrogenClass = E.Class;

if (E.Class < 0),                             % switch base is in standard position
  T     = E.NT1;
  E.NT1 = E.NT2;
  E.NT2 = T;
  T     = E.NT1.ID;
  E.NT1.ID = E.NT2.ID;
  E.NT2.ID = T;
  if E.NT1.Code ~= E.NT2.Code,
    E.Class = -E.Class;
  end
end

E.subplot = 4*(E.NT1.Code-1)+E.NT2.Code;      % which subplot to put this in

if ~isfield(E,'original'),
  E.original = 1;
end

ic = zIsostericSubgroups(E.NT1.Code,E.NT2.Code,abs(fix(E.Class)));

B(m) = E;                    % store this exemplar for isodisc calc
Lab{m} = [E.NT1.Base E.NT2.Base ' ' zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code) ic sprintf(' %5d ',E.Count) E.Source sprintf(' %5.1f',E.Class)];
Cat{m} = upper(zEdgeText(E.Class,Subcat,E.NT1.Code,E.NT2.Code));


return

zExemplarTable(1,0,0,1);
zExemplarTable(2,0,0,1);
zExemplarTable(3,0,0,1);
zExemplarTable(4,0,0,1);
zExemplarTable(5,0,0,1);
zExemplarTable(6,0,0,1);
zExemplarTable(7,0,0,1);
zExemplarTable(8,0,0,1);
zExemplarTable(9,0,0,1);
zExemplarTable(10,0,0,1);
zExemplarTable(11,0,0,1);
zExemplarTable(12,0,0,1);
zExemplarTable(13,0,0,1);
zExemplarTable([1:12],3.5,0,0);
zExemplarTable([1 5 4 7 10 11 2 6 3 8 9 12],3.5,0,0);

zExemplarTable([1 5],3.5,0);
zExemplarTable(1:6,3.5,0);
zExemplarTable(7:12,3.5,0);

% Neocles says:
% Separate the 2 groups: Group I: cWW, cWS, tWH, cHH, tHS, cSS
% Group II: tWW, tWS, cWH, tHH, cHS, tSS,

zExemplarTable(1,0,1,1);
zExemplarTable(2,0,1,1);
zExemplarTable(3,0,1,1);
zExemplarTable(4,0,1,1);
zExemplarTable(5,0,1,1);
zExemplarTable(6,0,1,1);
zExemplarTable(7,0,1,1);
zExemplarTable(8,0,1,1);
zExemplarTable(9,0,1,1);
zExemplarTable(10,0,1,1);
zExemplarTable(11,0,1,1);
zExemplarTable(12,0,1,1);
