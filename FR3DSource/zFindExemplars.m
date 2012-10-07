% zFindExemplars finds the best representative for each category of pairs.

Verbose = 1;              % Verbose = 1 tells it to show distance graphs
LMax    = 200;             % maximum number of pairs to consider in each class

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [6 7 13 14 15];
pcodes = [6 14 15 16];
pcodes = [6 7 9 11 13 14 15 16];    % pair codes to work on
pcodes = [1 5];
pcodes = [1 5 6 7 9 11 13 14 15 16];    % pair codes to work on

load('PairExemplars','Exemplar');

% clear Exemplar

CL = zClassLimits;

Pairs{1}  = 'AA';
Pairs{5}  = 'AC';
Pairs{6}  = 'CC';
Pairs{7}  = 'GC';
Pairs{9}  = 'AG';
Pairs{11} = 'GG';
Pairs{13} = 'AU';
Pairs{14} = 'CU';
Pairs{15} = 'GU';
Pairs{16} = 'UU';

% ------------------------------------------ Load non-redundant dataset

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData('NonRedundant_2008_02_21_list',2);   % load PDB data
  File = File(SIndex);
else
  [File,SIndex] = zAddNTData('NonRedundant_2008_02_21_list',2,File); % add PDB data if needed
  File = File(SIndex);
end                       

for f = 1:length(File),
  if ~isempty(File(f).Info.Resolution),
    Res(f) = File(f).Info.Resolution;
  else
    Res(f) = 10;
  end

  S(f) = length(File(f).NT);
end

[y,i] = sort(Res);
File = File(i);
S = S(i);

% ------------------------------------------ Load modeled basepairs

if ~exist('ModelFile'),                           % if no molecule data is loaded,
  [ModelFile,SIndex] = zAddNTData('Model_list',2);   % load PDB data
else
  [ModelFile,SIndex] = zAddNTData('Model_list',2,ModelFile); % add PDB data if needed
end                       

ModelFile = ModelFile(SIndex);

% ------------------------------------------ Load curated exemplars

if ~exist('CuratedFile'),                           % if no molecule data is loaded,
  [CuratedFile,SIndex] = zAddNTData('Curated_list',2);   % load PDB data
else
  [CuratedFile,SIndex] = zAddNTData('Curated_list',2,CuratedFile); % add PDB data if needed
end                       

CuratedFile = CuratedFile(SIndex);

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),            % run through all pair codes specified
 pc = pcodes(j);
 CLE = [CL(:,1,pc); [21 22 23]'];    % include stacking, identified differently
 CLE = CLE(find(CLE));               % leave out empty entries


 for row = 1:length(CLE),

  % specify criteria for selection of pairs ----------------------------------

  Param = [];
  Param(1,1) = CLE(row);
  Param(1,2) = pc;

  if abs(CLE(row)) < 20,
    Decimal = 1;
  else
    Decimal = 0;
  end

  % select pairs using selection criteria -----------------------------------

  List = zFindPairs(CuratedFile,Param,Decimal); % check curated basepairs

  if length(List) > 0,                    % this is a curated basepair

    fprintf('Using a curated pair for %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

%xDisplayCandidates(CuratedFile,List);

    f  = List(1,3);                       % file number
    Exemplar(row,pc).Filename   = CuratedFile(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = CuratedFile(f).NT(List(1,1));
    Exemplar(row,pc).NT2        = CuratedFile(f).NT(List(1,2));

    CList = zFindPairs(File,Param,Decimal);   % search the non-redundant list

    Exemplar(row,pc).Count      = length(CList(:,1));

  else

   List = zFindPairs(File,Param,Decimal);   % search the non-redundant list
   if length(List) > 0,               % instances of this category are found

    fprintf('Found %5d instances of %2s %4s class %5.1f ', length(List), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

    L = min(LMax,length(List(:,1))); % Limit the number of pairs to consider

% xDisplayCandidates(File,List);

    PD = zeros(L,L);               % Initialize the pair discrepancy
    ID = zeros(L,L);
    Ang = zeros(L,L);
    T1  = zeros(L,L);
    T2  = zeros(L,L);
    CP  = zeros(L,L);
    for k = 1:L,                   % Slow nested loop
      f1    = List(k,3);
      Model = List(k,[1 2]);
      NT1 = File(f1).NT(Model(1));
      NT2 = File(f1).NT(Model(2));
      for m = (k+1):L,
        f2    = List(m,3);
        Cand  = List(m,[1 2]);
        PD(k,m) = xDiscrepancy(File(f1),Model,File(f2),Cand);

        NT3 = File(f2).NT(Cand(1));
        NT4 = File(f2).NT(Cand(2));

        [d,ang,t1,t2,cp] = zIsoDiscrepancy(NT1,NT2,NT3,NT4);
        ID(k,m)  = d;
        Ang(k,m) = ang^2;
        T1(k,m)  = t1*t1';
        T2(k,m)  = t2*t2';
        CP(k,m)  = cp^2;
      end
    end

    PD = sqrt(PD)/2;            % finish discrepancy calculation

    bigm = max(max(PD));
    if bigm > 1,
      fprintf('%6.2f maximum pair discrepancy\n',bigm);
    else
      fprintf('\n');
    end

    PD = PD + PD';
    ID = ID + ID';

    rs = sum(PD);
    [y,i] = sort(rs);

    qs = sum(ID);
    [z,k] = sort(qs);

    if (length(List(:,1)) > 2) && (Verbose > 0),
      figure(1)
      clf
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
      FN = ['Instance_Discrepancies_' Pairs{pc} '_' zEdgeText(CLE(row),Decimal,pc)];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');

      drawnow

      figure(2)
      clf
      q = zClusterGraph(ID);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 10]);
      colorbar('location','eastoutside');
      hold on
      j = find(q == i(1));
      plot(j+0.5,j+0.5,'w*');
      m = find(q == k(1));
      plot(m+0.5,m+0.5,'wo');
      title(['IsoDiscrepancy between instances of ' Pairs{pc} ' ' zEdgeText(CLE(row),Decimal,pc)]);
      xlabel('Centroids marked (Geo discrep: white star, IDI: white circle)');
      FN = ['Instance_IsoDiscrepancies_' Pairs{pc} '_' zEdgeText(CLE(row),Decimal,pc)];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');


      drawnow
    end

    f = List(i(1),3);
    Exemplar(row,pc).Filename   = File(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = File(f).NT(List(i(1),1));
    Exemplar(row,pc).NT2        = File(f).NT(List(i(1),2));
    Exemplar(row,pc).Count      = length(List(:,1));

    else

      List = zFindPairs(ModelFile,Param,Decimal); % check modeled basepairs

      if length(List) > 0,

        fprintf('Using a model for       %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

        f = List(1,3);
        Exemplar(row,pc).Filename   = ModelFile(f).Filename;
        Exemplar(row,pc).Class      = CLE(row);
        Exemplar(row,pc).NT1        = ModelFile(f).NT(List(1,1));
        Exemplar(row,pc).NT2        = ModelFile(f).NT(List(1,2));
        Exemplar(row,pc).Count      = 0;
      else
        fprintf('No instances and no model for %2s %4s %6.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));
      end
    end
  end

  % add information to speed up the discrepancy calculation later

  [s,t] = size(Exemplar);
  if (row <= s) && (pc <= t),
   if ~isempty(Exemplar(row,pc).NT1),
    E = Exemplar(row,pc);
    Exemplar(row,pc).R            = E.NT2.Rot' * E.NT1.Rot;
    Exemplar(row,pc).T1           = (E.NT2.Center - E.NT1.Center) * E.NT1.Rot;
    Exemplar(row,pc).T2           = (E.NT1.Center - E.NT2.Center) * E.NT2.Rot;
    Exemplar(row,pc).AngleWeight  = [1 1];
    Exemplar(row,pc).LDiscCutoff  = Inf;
   end
  end 

  if LMax >= 500,
    save(['FR3DSource' filesep 'PairExemplars'],'Exemplar'); % Matlab version 7 only
    save PairExemplars_Version_6.mat Exemplar -V6 % for compatibility with older versions
  end

 end
end

zWriteExemplarPDB(1)
