% zExamplarTable(Cateogry) displays the best known
% representatives for interactions involving all pairs in
% interaction category(ies) Category
% Example:  zExemplarTable(1,3.5)
% Example:  zExemplarTable(1:12,3.5)

function [void] = zExemplarTable(Category,threshold)

if nargin == 1,
  threshold = 'default';
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
ViewParam.LabelBases= 8;                  % font size

% loop through computer classifications ----------------------

m = 1;

for ca = 1:length(Category),
 figure(fix(Category(ca)))
 clf
 plotted = zeros(16,1);

 for c1 = 1:4,
  for c2 = 1:4,
   for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar
    pc  = 4*(c2-1)+c1;                      % current paircode
    pc2 = 4*(c1-1)+c2;                      % current subplot number

    E = Exemplar(r,pc);                     % current exemplar

    if isempty(E.Filename),                 % no such exemplar
      E.Pair.Edge = 0;                          % avoid errors below
    end

    if (isempty(E.Filename) && any(Category(ca) == [1 2 7 8]) ...
       && any(pc == [2 3 4 8 10 12])),       % use symmetric case
      E     = Exemplar(r,pc2);
      E.NT1 = Exemplar(r,pc2).NT2;           % reverse nucleotides
      E.NT2 = Exemplar(r,pc2).NT1;
    elseif (E.Pair.Edge < 0) && ~(pc == 7 && any(Category(ca) == [1 2 7 8])),
      pc2   = pc;                            % display across diagonal
      T     = E.NT1;
      E.NT1 = E.NT2;
      E.NT2 = T;
      E.Pair.Edge = -E.Pair.Edge;
    end

    if ~isempty(E.Filename),


% modify the following line to include or exclude subcategories

%     if abs(E.Class) == Category(ca),          % don't include subcategories
 
     if abs(fix(E.Class)) == Category(ca),     % include subcategories

%       fprintf('%s%s %s %s%s %s Category %3.1f\n',E.NT1.Base,E.NT1.Number,E.Pair.EdgeText,E.NT2.Base,E.NT2.Number,E.Filename,E.Class);
%       zListPairData(E.Pair,1);

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
         xlab{pc2} = [xlab{pc2} ', ' num2str(E.Count)];        % append subcat count
       else
         ViewParam.LabelBases = 8;                  % font size
         xlab{pc2} = ['Count: ' num2str(E.Count)];
         Title{pc2} = [E.NT1.Base E.NT2.Base zEdgeText(E.Pair.Edge) ' ' strrep(E.Filename,'_','\_') ' '];
         Title{pc2} = [Title{pc2} E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
         CP = norm(E.NT1.Sugar(1,:) - E.NT2.Sugar(1,:));     % c1'-c1' dist
         Title{pc2} = [Title{pc2} ' ' num2str(CP)];
       end

       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;
       zDisplayNT(F,[1 2],ViewParam);
       zPlotHydrogenBonds(E.NT1,E.NT2,E.Class,E.NT1.Rot,E.NT1.Fit(1,:));

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

      if (E.Count >= 1) % && (fix(E.Class) == E.Class),
       B(m) = E;                                % store this exemplar
       Lab{m} = [E.NT1.Base E.NT2.Base zEdgeText(E.Pair.Edge) ' ' sprintf('%5.1f',E.Pair.Class) ' ' E.Filename ' ' sprintf('%4d',E.Count)];
       m = m + 1;

       if any(Category(ca) == [1 2 7 9 11 12]) && (E.NT1.Code == E.NT2.Code),
         B(m) = E;                              % store again
         B(m).NT1 = E.NT2;
         B(m).NT2 = E.NT1;
         B(m).Pair.Edge = -E.Pair.Edge;
         Lab{m} = [B(m).NT1.Base B(m).NT2.Base zEdgeText(B(m).Pair.Edge) ' ' sprintf('%5.1f',B(m).Pair.Class) ' ' B(m).Filename ' ' sprintf('%4d',B(m).Count) ];
         m = m + 1;
       end
      end

      end
    end
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
G = zeros((length(B)^2-length(B))/2,7);
j = 1;

for a = 1:length(B),
  for b = (a+1):length(B),

    R1 = B(a).NT2.Rot' * B(a).NT1.Rot;
    R2 = B(b).NT1.Rot' * B(b).NT2.Rot;

    ang = zAngleOfRotation(R1 * R2);

    t1 = (B(a).NT2.Sugar(1,:) - B(a).NT1.Sugar(1,:))*B(a).NT1.Rot - (B(b).NT2.Sugar(1,:) - B(b).NT1.Sugar(1,:))*B(b).NT1.Rot;
    t2 = (B(a).NT1.Sugar(1,:) - B(a).NT2.Sugar(1,:))*B(a).NT2.Rot - (B(b).NT1.Sugar(1,:) - B(b).NT2.Sugar(1,:))*B(b).NT2.Rot;

    cp = abs(norm(B(a).NT2.Sugar(1,:) - B(a).NT1.Sugar(1,:)) ...
           - norm(B(b).NT2.Sugar(1,:) - B(b).NT1.Sugar(1,:)));

    % calculate isodiscrepancy

    D(a,b) = sqrt((4*ang)^2 + (t1*t1' + t2*t2')/2 + (3*cp)^2);
    D(b,a) = D(a,b);

    G(j,:) = [6*ang sqrt(t1*t1') sqrt(t2*t2') 3*cp a b D(a,b)];
    j = j + 1;

  end
end

% ----------------------------------------- Cluster analysis and graph

Y = squareform(full(D));                       % convert to a vector

DD = full(D);

Z = linkage(Y,'average');                      % compute cluster tree
figure(fix(Category(ca))+13)
[H,T,p] = dendrogram(Z,0,'colorthreshold',threshold,'orientation','left','labels',Lab);

DDD = DD(p,p);
[s,t] = size(DDD);
if length(Category) == 1,
  fprintf('Table %d\n',Category(ca));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end
for i=1:s,
  fprintf('%24s ', Lab{p(i)});
  for j= 1:t,
    fprintf('%5.2f ', DDD(i,j));
  end
  fprintf('\n');
end
fprintf('\n');


h = gcf;
hh = gca;
if length(Lab) > 80,
  set(hh,'FontSize',3);
elseif length(Lab) > 40,
  set(hh,'FontSize',6);
else
  set(hh,'FontSize',8);
end
orient landscape
%print(h,'-depsc2',['Isostericity\ClusterIsoDisc' zEdgeText(Category)]);

% ----------------------------------------- Display nearby exemplars

% ------------ display in dendrogram order

[s,t] = size(DDD);
if length(Category) == 1,
  fprintf('Table %d\n',Category(1));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end

for i=1:s,
  fprintf('%24s ', Lab{p(i)});
  [d,h] = sort(DDD(i,:));
  a = 2;
  d = [d 99999];

  while (d(a) < 10),
    fprintf('%24s %5.2f ', Lab{p(h(a))}, d(a));
    a = a + 1;
  end

  fprintf('\n');
end
fprintf('\n');

% ------------ display in zClassLimits order

[s,t] = size(DD);
if length(Category) == 1,
  fprintf('Table %d\n',Category(1));
else
  fprintf('Tables ');
  for c = 1:length(Category),
    fprintf('%d ', Category(c));
  end
  fprintf('\n');
end

for i=1:s,
  fprintf('%24s ', Lab{i});
  [d,h] = sort(DD(i,:));
  a = 2;
  d = [d 99999];

  while (d(a) < 10),
    fprintf('%24s %5.2f ', Lab{h(a)}, d(a));
    a = a + 1;
  end
  fprintf('\n');
end
fprintf('\n');



return

% ----------------------------------------- Components of isodiscrepancy

Lab{1} = 'Angle';
Lab{2} = 'Vector diff 1';
Lab{3} = 'Vector diff 2';
Lab{4} = 'C1* difference';

g = 4;
figure(fix(Category(ca))+26)
clf
for i=2:g,
  for j=1:(i-1),
    subplot(g-1,g-1,(i-2)*(g-1)+j);
    plot(G(:,i),G(:,j),'.');
    xlabel(Lab{i});
    ylabel(Lab{j});
  end
end

end

% ----------------------------------------- Display basepairs by IsoDisc

fprintf('Comparison of basepairs, lowest IsoDiscrepancy first\n')

G = sortrows(G,7);
for i = 1:length(G(:,1)),
  a = G(i,5);
  b = G(i,6);
  figure(fix(Category(ca))+39)
  clf
  F.NT(1) = B(a).NT1;
  F.NT(2) = B(a).NT2;
  F.Filename = B(a).Filename;
  subplot(1,2,1)
  ViewParam.LineStyle = '-';
  zDisplayNT(F,[1 2],ViewParam);
  hold on

  subplot(1,2,2)
  hold on
  ViewParam.LineStyle = '--';
  zDisplayNT(F,[2 1],ViewParam);

  F.NT(1) = B(b).NT1;
  F.NT(2) = B(b).NT2;
  F.Filename = B(b).Filename;
  subplot(1,2,1)
  ViewParam.LineStyle = '--';
  zDisplayNT(F,[1 2],ViewParam);
  hold on
  view(2)
  subplot(1,2,2)
  rotate3d on
  ViewParam.LineStyle = '-';
  zDisplayNT(F,[2 1],ViewParam);
  hold on
  xlab = ['Count: ' num2str(B(b).Count)];
  Title = [B(b).NT1.Base B(b).NT2.Base zEdgeText(B(b).Pair.Edge) ' ' num2str(B(b).Pair.Class,3) ' ' strrep(B(b).Filename,'_','\_') ' '];
  Title = [Title B(b).NT1.Base B(b).NT1.Number '-' B(b).NT2.Base B(b).NT2.Number];
  CP = norm(B(b).NT1.Sugar(1,:) - B(b).NT2.Sugar(1,:));     % c1'-c1' dist
  Title = [Title ' ' num2str(CP)];
  xlabel(xlab);
  title(Title);
  view(2)

  subplot(1,2,1)
  rotate3d on
  xlab = ['Count: ' num2str(B(a).Count)];
  Title = [B(a).NT1.Base B(a).NT2.Base zEdgeText(B(a).Pair.Edge) ' ' num2str(B(a).Pair.Class,3) ' ' strrep(B(a).Filename,'_','\_') ' '];
  Title = [Title B(a).NT1.Base B(a).NT1.Number '-' B(a).NT2.Base B(a).NT2.Number];
  CP = norm(B(a).NT1.Sugar(1,:) - B(a).NT2.Sugar(1,:));     % c1'-c1' dist
  Title = [Title ' ' num2str(CP)];
  xlabel(xlab);
  title(Title);

  fprintf('IsoDiscrepancy %6.3f\n', G(i,7));
  fprintf('Angle contribution          %6.3f\n', G(i,1));
  fprintf('C1* location contribution 1 %6.3f\n', G(i,2));
  fprintf('C1* location contribution 2 %6.3f\n', G(i,3));
  fprintf('C1* distance contribution   %6.3f\n', G(i,4));
  fprintf('Press a key\n');

  pause
  
end
  
return

%zExemplarTable(1,3.5)
%zExemplarTable(2,3.5)
%zExemplarTable(3,3.5)
%zExemplarTable(4,3.5)
%zExemplarTable(5,3.5)
%zExemplarTable(6,3.5)
%zExemplarTable(7,3.5)
%zExemplarTable(8,3.5)
%zExemplarTable(9,3.5)
%zExemplarTable(10,3.5)
%zExemplarTable(11,3.5)
%zExemplarTable(12,3.5)
%zExemplarTable(13,3.5)
%zExemplarTable([1:13],3.5)
%zExemplarTable([1 5],3.5)
%zExemplarTable(1:6,3.5)
%zExemplarTable(7:12,3.5)

