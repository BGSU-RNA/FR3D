
figure(8)
clf
subplot(1,3,1)
cla
pp = OLO(PD);
q = zClusterGraph(PD,[],1,pp);
colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 1.4]);
colorbar('location','eastoutside');
hold on
j = find(q == i(1));
plot(j+0.5,j+0.5,'w*');
m = find(q == k(1));
plot(m+0.5,m+0.5,'wo');

title(['OLO Discrepancy for ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
%xlabel('Centroids (Geo discrep: star, IDI: circle)');
ylabel('Blue on diagonal means non coplanar');

subplot(1,3,2)
cla
pp = tpTSP(PD);
q = zClusterGraph(PD,[],1,pp);
colormap('default');
map = colormap;
map = map((end-8):-1:8,:);
colormap(map);
caxis([0 1.4]);
colorbar('location','eastoutside');
hold on
j = find(q == i(1));
plot(j+0.5,j+0.5,'w*');
m = find(q == k(1));
plot(m+0.5,m+0.5,'wo');

title(['tpTSP Discrepancy for ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
xlabel('Centroids (Geo discrep: star, IDI: circle)');
ylabel('Blue on diagonal means non coplanar');

subplot(1,3,3)
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
title(['IDI OLO Discrepancy for ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
%xlabel('Centroids (Geo discrep: star, IDI: circle)');
ylabel('Blue on diagonal means non coplanar');

%drawnow
FN = ['Instance_Discrepancies_' Pairs{pc} '_' zEdgeText(CLE(row),Decimal,pc)];
c = 8;
d = 2.47*c/6;
set(gcf, 'PaperPosition', [0 0 c d]); % 6 wide, 2.5 high

saveas(gcf,['Exemplars' filesep FN '.png'],'png');
