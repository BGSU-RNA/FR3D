% zFindBasepairExemplars finds the best representative for each category of pairs.

% File = zFindBasepairExemplars('1.53');          % the first time you run it in a Matlab session
% File = zFindBasepairExemplars('1.53',File);     % subsequent times, so you don't have to wait to load File again

% File = zFindBasepairExemplars('http://rna.bgsu.edu/rna3dhub/nrlist/download/current/4.0A/csv')

% function [File] = zFindBasepairExemplars(NRSetID,File)

tic

if ~(exist([pwd filesep 'FR3DSource']) == 7),        % if directory doesn't yet exist
  mkdir([pwd filesep 'FR3DSource']);
end

if ~(exist([pwd filesep 'Exemplars']) == 7),        % if directory doesn't yet exist
  mkdir([pwd filesep 'Exemplars']);
end

if nargin < 1 || length(NRSetID) == 0,
  NRSetID = 'http://rna.bgsu.edu/rna3dhub/nrlist/download/current/4.0A/csv';
end

if length(NRSetID) < 10,
  NRList = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/' NRSetID '/4.0A/csv'];
else
  NRList = NRSetID;
end

% Curated pairs are stored in PDBFiles and listed in Curated_list.pdb
% Every year or so, someone should set ViewCurated = 1 and run zFindExemplars
% It will display the curated exemplar in Figure 2 and show instances in
% the NR dataset in Figures 1 and 8.  Use xDisplayCandidates to view the
% top candidates.  If they are as satisfactory as the curated exemplar,
% remove the curated exemplar from Curated_list.pdb.
% CLZ did this on 2011-07-26.
% CLZ did this on 2013-10-19.
% 2013-10-19 removed Curated_Model_cWS_UU_Exemplar
% 2013-10-19 removed Curated_1J5E_cSs_GG_Exemplar
% 2013-10-19 removed Curated_2J01_cHS_GA_Exemplar

Verbose = 1;                     % Verbose = 1 tells it to show distance graphs
ViewCurated = 0;                 % ViewCurated = 1 stops to review each curated pair, but *only* curated pairs
ViewInstancesOfModeledPairs = 1; % Stop when an instance of a formerly modeled pair is found
LMax    = 500;                   % maximum number of pairs to consider in each class
UseCrystalSymmetries = 1;        % use .pdb1 files?
VerifyCutoffs = 1;               % stop after each motif group to interactively view instances

% Pair codes:  1-AA 5-AC 6-CC 7-GC 9-AG 11-GG 13-AU 14-CU 15-GU 16-UU

pcodes = [9 11 13 14 15 16 1 5 6 7 ];    % pair codes to work on

load(['PairExemplars'],'Exemplar');       % load previous exemplars

% ------------------------------------- % consistency check
for c1 = 1:4,
  for c2 = 1:4,
     pc  = 4*(c2-1)+c1;                       % current paircode
     for r = 1:length(Exemplar(:,1)),         % loop through rows of Exemplar
        E = Exemplar(r,pc);                    % current exemplar
        if ~isempty(E.NT1),                    % non-empty entry of Exemplar
           Epc = 4*(E.NT2.Code-1)+E.NT1.Code;     % paircode of exemplar itself; some things are wrong!
           if Epc ~= pc,
             fprintf('Paircode disagreement, %d versus %d in row %d\n', pc, Epc, r);
           end
       end
      end
   end
end

% load('PairExemplars_Old','Exemplar');   % load last established

OldExemplar = Exemplar;                 % for when no instances are found

Exemplar = [];                          % start fresh

CL = zClassLimits;                      % load basepair classification limits

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
  [File,SIndex] = zAddNTData(NRList,0,[],1);   % load PDB data
  File = File(SIndex);
else
  [File,SIndex] = zAddNTData(NRList,0,File,1); % add PDB data if needed
  File = File(SIndex);
end

% ------------------------------------------ Omit .pdb1, .pdb2 files?

if UseCrystalSymmetries == 0,
  OK = ones(1,length(File));
  for f = 1:length(File),
    if File(f).PDBFilename(end) == '1',
      OK(f) = 0;
    end
  end
  j = find(OK);
  fprintf('Omitting %d structures that use crystal symmetries\n',length(File) - length(j));
  File = File(j);
end

% ------------------------------------------- Gather resolution information and sort files high res to low res

for f = 1:length(File),
  if ~isempty(File(f).Info.Resolution),
    Res(f) = File(f).Info.Resolution;
  else
    Res(f) = 4;                                % unknown, but we use a 4A list, so it should not be worse than 4A
  end

  S(f) = length(File(f).NT);                   % size of each file
end

[y,i] = sort(Res);                             % sort files by resolution
File = File(i);
Res = Res(i);
S = S(i);

% ------------------------------------------ Load modeled basepairs
% Note:  If the file names in Model_list.pdb end with .pdb, it will load the pdb file, which is slow
%        If they do not end with .pdb, it will read .mat files, if they exist, but if they don't, it will look for .cif files, which it won't find

if ~exist('ModelFile'),                        % if no molecule data is loaded,
  [ModelFile,SIndex] = zAddNTData('Model_list',0,[],1);   % load PDB data
else
  [ModelFile,SIndex] = zAddNTData('Model_list',0,ModelFile,1);
end

for i = 1:length(ModelFile),
  ModelFile(i).Info.ExpTechnique = 'Modeled';
  zSaveNTData(ModelFile(i));
end

ModelFile = ModelFile(SIndex);

% ------------------------------------------ Load curated exemplars

if ~exist('CuratedFile'),                      % if no molecule data is loaded,
  [CuratedFile,SIndex] = zAddNTData('Curated_list',0,[],1);   % load PDB data
else
  [CuratedFile,SIndex] = zAddNTData('Curated_list',0,CuratedFile,1);
end

for i = 1:length(ModelFile),
  CuratedFile(i).Info.ExpTechnique = 'Modeled';
  zSaveNTData(CuratedFile(i));
end

CuratedFile = CuratedFile(SIndex);

% ------------------------------------------ Make a directory for exemplars

if ~(exist('Exemplars') == 7),        % if directory doesn't yet exist
  mkdir('Exemplars');
end

% loop through paircodes and computer classifications ----------------------

for j = 1:length(pcodes),            % run through all pair codes specified
 pc = pcodes(j);                     % current paircode
fprintf('Current paircode is %2d\n',pc);
 CLE = CL(:,1,pc);                   % class limits for this paircode
 CLE = CLE(find(CLE));               % leave out empty entries
 CLE = CLE(find(abs(CLE) < 20));     % leave out stacking

 for row = 1:length(CLE),            % run through classes for this paircode

  % specify criteria for selection of pairs ----------------------------------

  Param = [];
  Param(1,1) = CLE(row);
  Param(1,2) = pc;

  if abs(CLE(row)) < 20,
    Decimal = 1;                          % treat subcategories separately
  else
    Decimal = 0;                          % group subcategories together
  end

  % select pairs using selection criteria -----------------------------------

  List = zFindPairs(CuratedFile,Param,Decimal); % check curated basepairs

  if length(List) > 0,                    % this is a curated basepair

    fprintf('Using a curated pair for %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

    f  = List(1,3);                       % file number
    Exemplar(row,pc).Filename   = CuratedFile(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).NT1        = CuratedFile(f).NT(List(1,1));
    Exemplar(row,pc).NT2        = CuratedFile(f).NT(List(1,2));

    Exemplar(row,pc).Status     = 'Curated';

    if strcmp(lower(CuratedFile(f).Filename(9:12)),'mode'),
      Exemplar(row,pc).Source = 'Model';
      Exemplar(row,pc).Resolution = NaN;
      Exemplar(row,pc).NT1ID      = '';
      Exemplar(row,pc).NT2ID      = '';
      Exemplar(row,pc).NT1.Number = '1';
      Exemplar(row,pc).NT2.Number = '2';
      Exemplar(row,pc).NT1.Chain  = 'A';
      Exemplar(row,pc).NT2.Chain  = 'A';
    else
      FN = upper(CuratedFile(f).Filename(9:12));
      Exemplar(row,pc).Source = FN;    % extract PDB ID
      NewFile = zAddNTData(FN);
      Exemplar(row,pc).Resolution = NewFile.Info.Resolution;
      i1 = zIndexLookup(NewFile,Exemplar(row,pc).NT1.Number,Exemplar(row,pc).NT1.Chain);
      Exemplar(row,pc).NT1ID      = aGetNTId_NDB(NewFile,i1);
      i2 = zIndexLookup(NewFile,Exemplar(row,pc).NT2.Number,Exemplar(row,pc).NT2.Chain);
      Exemplar(row,pc).NT2ID      = aGetNTId_NDB(NewFile,i2);
%fprintf('http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=%s,%s\n',Exemplar(row,pc).NT1ID,Exemplar(row,pc).NT2ID);
    end

    List = zFindPairs(File,Param,Decimal);   % search the non-redundant list

    Exemplar(row,pc).Count = 0;

    if length(List) > 0,               % instances of this category are found in 3D structures now
      Exemplar(row,pc).Count = length(List(:,1));

      fprintf('FYI, found %5d instances of %2s %4s class %5.1f\n', length(List), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));

      [PD,ID,i] = zOrderPairs(File,List,LMax,Verbose);

      rs = sum(PD);
      gg = Res(List(i(1:min(end,80)),3))';             % resolution of the structure
      hh = rs(i(1:min(end,80)))';                      % row sum of discrepancies
      kk = 1+((1:length(rs))/length(rs))';             % 1 + rank among instances

      [y,w] = sort(gg.*hh.*kk);

      if 0 > 1,
        disp('Top-ranking candidates, last column is the overall criterion')
        [gg(1:min(end,20)) hh(1:min(end,20)) kk(1:min(end,20)) gg(1:min(end,20)).*hh(1:min(end,20)).*kk(1:min(end,20))]
        disp('Chosen candidate')
        [gg(w(1)) hh(w(1)) gg(w(1))*hh(w(1))*kk(w(1))]
      end

      i = i(w);

      NT1 = Exemplar(row,pc).NT1;
      NT2 = Exemplar(row,pc).NT2;
      f   = List(i(1),3);
      NT3 = File(f).NT(List(i(1),1));
      NT4 = File(f).NT(List(i(1),2));
      [File(f).Filename ' ' num2str(File(f).Info.Resolution) ' ' NT3.Base NT3.Number ' ' NT4.Base NT4.Number]
      d   = zIsoDiscrepancy(NT1,NT2,NT3,NT4);
      if Verbose > 0,
        fprintf('    IsoDiscrepancy between curated and best from search is %7.4f\n',d);
      end

      if ViewCurated > 0,
        if (length(List(:,1)) > 2) && (Verbose > 0),
          zFindExemplarsPlot
        end
        close all
        figure(2)
        clf
        VP.Sugar = 1;
        VP.AtOrigin = 1;
        VP.LabelBases = 10;
        clear FFF
        FFF.NT(1) = Exemplar(row,pc).NT1;
        FFF.NT(2) = Exemplar(row,pc).NT2;
        zDisplayNT(FFF,1:2,VP);
        view(2)
        xlabel(['Curated pair from ' strrep(Exemplar(row,pc).Filename,'_','-')]);
        axis equal
        rotate3d on

        figure(1)
        clf
        xDisplayCandidates(File,List(i,:));
        rotate3d on
      end

    end

  elseif ViewCurated == 0,                  % not stopping to look at curated basepairs

   List = zFindPairs(File,Param,Decimal);   % search the non-redundant list

   if length(List) > 0,               % instances of this category are found

    fprintf('Found %5d instances of %2s %4s class %5.1f\n', length(List), Pairs{pc}, zEdgeText(CLE(row),Decimal,pc), CLE(row));
    if length(List) > LMax,
      fprintf('Restricting to the first %d of them from the highest resolution structures\n',LMax);
    end

    [PD,ID,i,k] = zOrderPairs(File,List,LMax,Verbose); % order by centrality

    rs = sum(PD);
    gg = Res(List(i(1:min(end,LMax)),3))';
    hh = rs(i(1:min(end,LMax)))';
    kk = 1+(((1:length(hh)))/length(rs))';          % rank among instances

    % the chosen instance minimizes the product of resolution, mean distance to other instances, and rank among instances
    % note that coplanar basepairs have a much lower value of PD and rs than non-coplanar; they are in different categories

    [y,w] = sort(gg.*hh.*kk);

    if 0 > 1,
      disp('Resolution, centrality, and rank of the most central instances:');
      [gg(1:min(end,20)) hh(1:min(end,20)) kk(1:min(end,20)) gg(1:min(end,20)).*hh(1:min(end,20)).*kk(1:min(end,20))]
      [gg(w(1)) hh(w(1)) gg(w(1))*hh(w(1))*kk(w(1))]
    end

    i = i(w);

    if (length(List(:,1)) > 2) && (Verbose > 0),
      zFindExemplarsPlot                          % display mutual discreps
    end

    f = List(i(1),3);
    Exemplar(row,pc).Filename   = File(f).Filename;
    Exemplar(row,pc).Class      = CLE(row);
    Exemplar(row,pc).Status     = 'Current_NR_Set';
    Exemplar(row,pc).Source     = File(f).Filename;
    Exemplar(row,pc).NT1        = File(f).NT(List(i(1),1));
    Exemplar(row,pc).NT2        = File(f).NT(List(i(1),2));
    Exemplar(row,pc).NT1ID      = aGetNTId_NDB(File(f),List(i(1),1));
    Exemplar(row,pc).NT2ID      = aGetNTId_NDB(File(f),List(i(1),2));
%fprintf('http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=%s,%s\n',Exemplar(row,pc).NT1ID,Exemplar(row,pc).NT2ID);
    Exemplar(row,pc).Count      = length(List(:,1));
    Exemplar(row,pc).Resolution = File(f).Info.Resolution;
    Exemplar(row,pc).PDBFilename= File(f).PDBFilename;

    MList = zFindPairs(ModelFile,Param,Decimal); % check modeled basepairs

    if length(MList) > 0,
      fprintf('Found an example of a modeled basepair!\n')
      if ViewInstancesOfModeledPairs > 0,
        figure(2)
        clf
        VP.Sugar = 1;
        VP.AtOrigin = 1;
        VP.LabelBases = 10;
        clear FFF
        FFF.NT(1) = Exemplar(row,pc).NT1;
        FFF.NT(2) = Exemplar(row,pc).NT2;
        zDisplayNT(FFF,1:2,VP);
        view(2)
        xlabel(['Curated pair from ' strrep(Exemplar(row,pc).Filename,'_','-')]);
        axis equal
        rotate3d on

        figure(1)
        clf
        xDisplayCandidates(File,List(i,:));
        rotate3d on
      end
    end

    if VerifyCutoffs > 0,
      figure(1)
      close
      figure(1)
      clf
      xDisplayCandidates(File,List(i,:));
      rotate3d on
    end

   else                                           % no instance found

      List = zFindPairs(ModelFile,Param,Decimal); % check modeled basepairs

      if length(List) > 0,                        % model was found

        fprintf('Using a model for       %2s %4s class %5.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

        f = List(1,3);
        Exemplar(row,pc).Filename   = ModelFile(f).Filename;
        Exemplar(row,pc).Class      = CLE(row);
        Exemplar(row,pc).Status     = 'Model';
        Exemplar(row,pc).Source     = 'Model';
        Exemplar(row,pc).NT1        = ModelFile(f).NT(List(1,1));
        Exemplar(row,pc).NT2        = ModelFile(f).NT(List(1,2));
        Exemplar(row,pc).NT1.Number = '1';
        Exemplar(row,pc).NT2.Number = '2';
        Exemplar(row,pc).NT1.Chain  = 'A';
        Exemplar(row,pc).NT2.Chain  = 'A';

        Exemplar(row,pc).NT1ID      = '';
        Exemplar(row,pc).NT2ID      = '';
        Exemplar(row,pc).Count      = 0;
        Exemplar(row,pc).Resolution = NaN;
        Exemplar(row,pc).PDBFilename= ModelFile(f).Filename;

      else
        fprintf('No instances and no model for %2s %4s %6.1f\n', Pairs{pc}, zEdgeText(CLE(row),1), CLE(row));

        Code1 = mod(pc-1,4)+1;
        Code2 = floor((pc-1)/4)+1;

        if CLE(row) == fix(CLE(row)),

          [N1,N2,E] = zGetExemplar(CLE(row),Code1,Code2,OldExemplar);

          opc = 4*(N2.Code-1) + N1.Code;        % verify old paircode!!

          if ~isempty(E.Filename) && opc == pc,
            fprintf('Using previous exemplar for this class\n');
            Exemplar(row,pc).Filename   = E.Filename;
            Exemplar(row,pc).Class      = E.Class;
            if length(E.Filename) == 4,
              NewFile = zAddNTData(E.Filename);
              Exemplar(row,pc).Status     = 'Previous_NR_Set';
              Exemplar(row,pc).Source     = E.Filename;
              Exemplar(row,pc).Resolution = NewFile.Info.Resolution;
              i1 = zIndexLookup(NewFile,E.NT1.Number,E.NT1.Chain);
              Exemplar(row,pc).NT1ID      = aGetNTId_NDB(NewFile,i1);
              i2 = zIndexLookup(NewFile,E.NT2.Number,E.NT2.Chain);
              Exemplar(row,pc).NT2ID      = aGetNTId_NDB(NewFile,i2);
%  fprintf('http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=%s,%s\n',Exemplar(row,pc).NT1ID,Exemplar(row,pc).NT2ID);
            else

  fprintf('This case should not happen ********************');

              Exemplar(row,pc).Status     = 'Old_Model';
              Exemplar(row,pc).Resolution = NaN;
              Exemplar(row,pc).NT1ID      = '';
              Exemplar(row,pc).NT2ID      = '';
            end
            Exemplar(row,pc).NT1        = E.NT1;
            Exemplar(row,pc).NT2        = E.NT2;

            Exemplar(row,pc).Count      = E.Count;
            Exemplar(row,pc).PDBFilename= E.PDBFilename;

          else
            fprintf('***** No previous exemplar found, may be time to make a new one.\n')
          end
        else
          fprintf('This is a basepair subcategory and may no longer be actively maintained, so no exemplar will be generated\n');
        end
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

  if LMax >= 500 && ViewCurated == 0,
    save([pwd filesep 'FR3DSource' filesep 'PairExemplars'],'Exemplar'); % Matlab version 7 only
    save PairExemplars_Version_6.mat Exemplar -V6 % for compatibility with older versions
  end

 end
end

fprintf('Note:  Copy PairExemplars.mat from local FR3DSource to Github FR3DSource\n')

% ----------------------------------------- Fix model fields

[s,t] = size(Exemplar);

clear NewExemplar

for i = 1:s,
  for j = 1:t,
    if ~isempty(Exemplar(i,j).Filename),
      E = Exemplar(i,j);
      if isfield(E.NT1,'ModelNum'),
        E.NT1.Model = E.NT1.ModelNum;
        E.NT2.Model = E.NT2.ModelNum;
        E.NT1 = rmfield(E.NT1,'ModelNum');
        E.NT2 = rmfield(E.NT2,'ModelNum');
      end
      if ~isfield(E.NT1,'Model'),
        E.NT1.Model = 1;
        E.NT2.Model = 1;
      end
      if ~isfield(E.NT1,'Unit'),
        E.NT1.Unit = E.NT1.Base;
        E.NT2.Unit = E.NT2.Base;
      end
      E.NT1 = orderfields(E.NT1);
      E.NT2 = orderfields(E.NT2);
      E = orderfields(E);
      NewExemplar(i,j) = E;
    end
  end
end

Exemplar = NewExemplar;

if LMax >= 500 && ViewCurated == 0,
  save([pwd filesep 'FR3DSource' filesep 'PairExemplars'],'Exemplar'); % Matlab version 7 only
  save PairExemplars_Version_6.mat Exemplar -V6 % for compatibility with older versions
end

% ------------------------------------------- Output results

zWriteExemplarPDB(Exemplar,1)

zExemplarIDICalculation                    % calculate and store ExemplarIDI
zExemplarFrequencyCalculation              % calculate and store ExemplarFreq

figure(1)
clf
figure(2)
clf

clear FN

% -------- Create isostericity tables without subcategories
[T{1},U,FN] = zExemplarTable(1,0,0,2);
AllFN = [FN];
[T{2},U,FN] = zExemplarTable(2,0,0,2);
AllFN = [AllFN FN];
[T{3},U,FN] = zExemplarTable(3,0,0,2);
AllFN = [AllFN FN];
[T{4},U,FN] = zExemplarTable(4,0,0,2);
AllFN = [AllFN FN];
[T{5},U,FN] = zExemplarTable(5,0,0,2);
AllFN = [AllFN FN];
[T{6},U,FN] = zExemplarTable(6,0,0,2);
AllFN = [AllFN FN];
[T{7},U,FN] = zExemplarTable(7,0,0,2);
AllFN = [AllFN FN];
[T{8},U,FN] = zExemplarTable(8,0,0,2);
AllFN = [AllFN FN];
[T{9},U,FN] = zExemplarTable(9,0,0,2);
AllFN = [AllFN FN];
[T{10},U,FN] = zExemplarTable(10,0,0,2);
AllFN = [AllFN FN];
[T{11},U,FN] = zExemplarTable(11,0,0,2);
AllFN = [AllFN FN];
[T{12},U,FN] = zExemplarTable(12,0,0,2);
AllFN = [AllFN FN];

% append_pdfs(['Isostericity' filesep 'All_basepairs.pdf'],AllFN{:});  % Make a nice PDF file that Craig likes; also uncomment squeeze_axes in zExemplarTable

zExemplarTable(13,0,0,1);
zExemplarTable(14,0,0,1);
zExemplarTable(15,0,0,1);

zExemplarFrequencyTables(T)

toc

return

zExemplarTablesExcel

return

zExemplarTable(1,0,1,2);
zExemplarTable(2,0,1,2);
zExemplarTable(3,0,1,2);
zExemplarTable(4,0,1,2);
zExemplarTable(5,0,1,2);
zExemplarTable(6,0,1,2);
zExemplarTable(7,0,1,2);
zExemplarTable(8,0,1,2);
zExemplarTable(9,0,1,2);
zExemplarTable(10,0,1,2);
zExemplarTable(11,0,1,2);
zExemplarTable(12,0,1,2);

