% zAmalTripleStudy is a script that can be run from xDisplayCandidates.m

    N1 = File(f).NT(Indices(1));
    N2 = File(f).NT(Indices(2));
    N3 = File(f).NT(Indices(3));

    if isfield(Search,'AvgDisc'),   
      ylabel(['Plot ',int2str(n),' of ',int2str(s),'   Average discrepancy from others ', num2str(Search.AvgDisc(n))]);
    elseif Query.Geometric > 0,
      ylabel(['Plot ',int2str(n),' of ',int2str(s),'   Discrepancy ',...
          num2str(Search.Discrepancy(n))]);
    else
      ylabel(['Plot ',int2str(n),' of ',int2str(s)]);
    end

    BP12 = File(f).Edge(Indices(1),Indices(2));
    BP23 = File(f).Edge(Indices(2),Indices(3));
    BP31 = File(f).Edge(Indices(3),Indices(1));

    Nam = [N1.Base N2.Base N3.Base '_' zEdgeText(BP12,1) '_' zEdgeText(BP23,1) '_' zEdgeText(BP31,1) '_Model_'];

    if BP12 > 100,  BP12 = BP12 - 100; end
    if BP12 < -100, BP12 = BP12 + 100; end
    if BP23 > 100,  BP23 = BP23 - 100; end
    if BP23 < -100, BP23 = BP23 + 100; end
    if BP31 > 100,  BP31 = BP31 - 100; end
    if BP31 < -100, BP31 = BP31 + 100; end

    if exist('comparetoregular.txt') > 0,
      BP12 = fix(BP12);
      BP23 = fix(BP23);
      BP31 = fix(BP31);
      Nam = [Nam zEdgeText(BP12,1) '_' zEdgeText(BP23,1) '_' zEdgeText(BP31,1) '_Model_'];
    end

    if BP12 ~= 0 && BP23 ~= 0 && BP31 ~= 0,
      FFF = zMakeTriple(BP12,BP23,N1.Code,N2.Code,N3.Code);
      if ~isempty(FFF),
        hold on
        VPP.AtOrigin = 2;
        VPP.LabelBases = 0;
        VPP.Sugar = 0;
        VPP.LineStyle = '-.';
        VPP.Title = 0;
        zDisplayNT(FFF,1:3,VPP)
    ytext = 'Model interactions ';
    ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
    ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
%    ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];
    xlabel(ytext);
        view(2)

  saveas(gcf,['Triples' filesep Nam '.fig'],'fig');
  saveas(gcf,['Triples' filesep Nam '.png'],'png');

        fprintf('12-23\n');
        fprintf('Geometric discrepancy from model: %7.4f\n',xDiscrepancy(File(f),Indices(1:3),FFF,1:3));
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3));

      end

      figure(34)
      clf
      FFF = zMakeTriple(BP31,BP12,N3.Code,N1.Code,N2.Code);
      if ~isempty(FFF),
        hold on
        FFF.NT = FFF.NT([2 3 1]);
        VPP.AtOrigin = 2;
        VPP.LabelBases = 0;
        VPP.Sugar = 0;
        VPP.LineStyle = '-.';
        VPP.Title = 0;
        zDisplayNT(FFF,1:3,VPP)
    ytext = 'Model interactions ';
    ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
%    ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
    ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];
    xlabel(ytext);
        view(2)
  axis vis3d
  axis equal

  saveas(gcf,['Triples' filesep Nam 'alternative.fig'],'fig');
  saveas(gcf,['Triples' filesep Nam 'alternative.png'],'png');


        fprintf('12-13\n');
        fprintf('Geometric discrepancy from model: %7.4f\n',xDiscrepancy(File(f),Indices(1:3),FFF,1:3));
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3));
      end

    elseif BP12 ~= 0 && BP23 ~= 0,
      FFF = zMakeTriple(BP12,BP23,N1.Code,N2.Code,N3.Code);
      if ~isempty(FFF),
        hold on
        VPP.AtOrigin = 2;
        VPP.LabelBases = 0;
        VPP.Sugar = 0;
        VPP.LineStyle = '-.';
        VPP.Title = 0;
        zDisplayNT(FFF,1:3,VPP)
    ytext = 'Model interactions ';
    ytext = [ytext ' 12:' zEdgeText(BP12,1) ' '];
    ytext = [ytext ' 23:' zEdgeText(BP23,1) ' '];
%    ytext = [ytext ' 13:' zEdgeText(BP31,1) ' '];
    xlabel(ytext);
        view(2)

  saveas(gcf,['Triples' filesep Nam '.fig'],'fig');
  saveas(gcf,['Triples' filesep Nam '.png'],'png');

        fprintf('12-23\n');
        fprintf('Geometric discrepancy from model: %7.4f\n',xDiscrepancy(File(f),Indices(1:3),FFF,1:3));
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3));

      end

    elseif BP12 ~= 0 && BP31 ~= 0,
      FFF = zMakeTriple(BP31,BP12,N3.Code,N1.Code,N2.Code);
      if ~isempty(FFF),
        hold on
        FFF.NT = FFF.NT([2 3 1]);
        VPP.AtOrigin = 2;
        VPP.LabelBases = 0;
        VPP.Sugar = 0;
        VPP.LineStyle = '-.';
        VPP.Title = 0;
        zDisplayNT(FFF,1:3,VPP)
    ytext = 'Model interactions ';
    ytext = [ytext ' 12:' zEdgeText(File(f).Edge(Indices(1),Indices(2))) ' '];
%    ytext = [ytext ' 23:' zEdgeText(File(f).Edge(Indices(2),Indices(3))) ' '];
    ytext = [ytext ' 13:' zEdgeText(File(f).Edge(Indices(1),Indices(3))) ' '];
    xlabel(ytext);
        view(2)

  saveas(gcf,['Triples' filesep Nam '.fig'],'fig');
  saveas(gcf,['Triples' filesep Nam '.png'],'png');

        fprintf('12-13\n');
        fprintf('Geometric discrepancy from model: %7.4f\n',xDiscrepancy(File(f),Indices(1:3),FFF,1:3));
        fprintf('Triple IsoDiscrepancy from model: %7.4f\n',xDiscrepancyForTriples(File(f),Indices(1:3),FFF,1:3));
      end
    end
