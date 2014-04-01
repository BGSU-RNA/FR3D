
      figure(8)
      clf
      subplot(1,2,1)
      cla
      q = zClusterGraph(PD);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 1]);
      colorbar('location','eastoutside');
      hold on
      j = find(q == i(1));
      plot(j+0.5,j+0.5,'w*');
      m = find(q == k(1));
      plot(m+0.5,m+0.5,'wo');

      title(['Discrepancy between instances of ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
      xlabel('Centroids marked (Geo discrep: white star, IDI: white circle)');
      ylabel('Blue on diagonal means non coplanar');
      FN = ['Instance_Discrepancies_' Pairs{pc} '_' zEdgeText(CLE(row),Decimal,pc)];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');

      drawnow

      subplot(1,2,2)
      cla
      q = zClusterGraph(ID);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 5]);
      colorbar('location','eastoutside');
      hold on
      j = find(q == i(1));
      plot(j+0.5,j+0.5,'w*');
      m = find(q == k(1));
      plot(m+0.5,m+0.5,'wo');
      title(['IsoDiscrepancy between instances of ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
      xlabel('Centroids marked (Geo discrep: white star, IDI: white circle)');
      ylabel('Blue on diagonal means non coplanar');
      FN = ['Instance_IsoDiscrepancies_' Pairs{pc} '_' zEdgeText(CLE(row),Decimal,pc)];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');


      drawnow
