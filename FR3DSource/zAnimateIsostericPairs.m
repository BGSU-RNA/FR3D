% zAnimateIsostericPairs makes bmp's that can be combined together into animated gif's of various basepairs, to show isostericity between them


load PairExemplars

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

clear Pairs

Pairs{1} = [13 1  1]; % paircode class display interchanged or not
Pairs{2} = [ 7 1  1];
Pairs{3} = [ 7 1 -1];
Pairs{4} = [13 1 -1];
ViewParam.AnimationFilename = 'cWW';

clear Pairs

Pairs{1} = [ 9 10  1];
Pairs{2} = [ 1 10  1];
Pairs{3} = [ 5 10  1];
Pairs{4} = [13 10  1];
Pairs{5} = [ 5 -10  -1];
Pairs{6} = [ 6 10  1];
Pairs{7} = [14 10  1];
ViewParam.AnimationFilename = 'tHS';

clear Pairs

Pairs{1} = [ 1 5  1];
Pairs{2} = [ 5 5  1];
Pairs{3} = [ 9 5  1];
Pairs{4} = [13 5  1];
ViewParam.AnimationFilename = 'cWS';

clear Pairs

Pairs{1} = [13 -4  -1];
Pairs{2} = [ 6  4   1];
ViewParam.AnimationFilename = 'tWH';

clear Pairs

Pairs{1} = [ 1 12  1];
Pairs{2} = [ 5 12  1];
Pairs{3} = [ 9 12  1];
Pairs{4} = [13 12  1];
ViewParam.AnimationFilename = 'tSs';

ViewParam.AtOrigin = 1;
ViewParam.Sugar    = 1;
ViewParam.LabelBases = 0;
ViewParam.Grid     = 0;
ViewParam.Animate  = 0;

for j = 1:length(Pairs),
   pc = Pairs{j}(1);

   for row = 1:length(Exemplar(:,pc)),

     E = Exemplar(row,pc);

     if ~isempty(E),
      if E.Class == Pairs{j}(2),

       Pair = zClassifyPair(E.NT1,E.NT2);
       fprintf('%s%s %s %s%s %s Category %3.1f\n',E.NT1.Base,E.NT1.Number,Pair.EdgeText,E.NT2.Base,E.NT2.Number,E.Filename,E.Class);
       zListPairData(Pair,1);

% display the exemplar pair ---------------------------------------------

       figure(1)
       clf

       F.NT(1) = E.NT1;
       F.NT(2) = E.NT2;
       F.Filename = E.Filename;

       if Pairs{j}(3) > 0,
         zDisplayNT(F,[1 2],ViewParam);
         zPlotHydrogenBonds(E.NT1,E.NT2,E.Class,E.NT1.Rot,E.NT1.Fit(1,:));
       else
         zDisplayNT(F,[2 1],ViewParam);
         zPlotHydrogenBonds(E.NT2,E.NT1,E.Class,E.NT2.Rot,E.NT2.Fit(1,:));
       end

       view(2)
       axis equal

       switch fix(E.Class),
         case {1, 2},    axis([-7 12 -4 11 -5 5]);
         case {4, -4},   axis([-7 10 -4 12 -5 5]);
         case {5, -5},   axis([-7 12 -3 13 -5 5]);
         case {10, -10}, axis([-8  6 -4 14 -5 5]);
         case {12, -12}, axis([-6  13 -5 9 -5 5]);
         case 15,        axis([-5 3 -3 5 -5 5]);
         case 16,        axis([-4 4 -3 5 -3 5]);
         otherwise,      axis([-6 10 -6 10 -6 10]);
       end

       Title = [E.NT1.Base E.NT2.Base ' ' num2str(E.Class) ' ' strrep(E.Filename,'_','\_') ' '];
       Title = [Title E.NT1.Base E.NT1.Number '-' E.NT2.Base E.NT2.Number];
       title(Title);

    set(gca,'TickLength',[0 0]);
    Image = getframe;
    Image.cdata = Image.cdata(1:(end-1),2:end,:);
    P = frame2im(Image);
    filename = [ViewParam.AnimationFilename num2str(j+10) '.bmp'];
    imwrite(P,filename, 'bmp');

end
end
end

%pause

end
