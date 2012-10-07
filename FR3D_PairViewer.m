% PairViewer(File,Param,ViewParam)
% loads one or more pdb files, allows the user to specify
% selection criteria, then displays pair information in various formats

function [File] = PairViewer(File,Param,ViewParam)

if nargin == 0,                  % If files haven't already been loaded
  Filenames = {'1s72'};
  File = zAddNTData(Filenames,0);
end

if exist('ViewParam'),
  if ~isfield(ViewParam,'FigNum'),
    ViewParam.FigNum = 2;
  end
else
  ViewParam.FigNum = 2;
end

while ViewParam.FigNum > 0,
  if nargin <= 1,
    if exist('Param'),
      Param = zEnterPairSelection(Param);
    else
      Param = zEnterPairSelection([]);
    end
  end

  SP = zSelectPairs(File,Param);

  if length(SP) > 0,
    if nargin < 3,
      ViewParam.Mode = 1;
    end
    while ViewParam.Mode(1) > 0,
      if nargin < 3,
        ViewParam = zEnterViewMode(Param,ViewParam);
      end

      if any(ViewParam.Mode == 1) || any(ViewParam.Mode == 6),
        SP = zColorPairs(File,SP,Param,ViewParam);
      end

      if ViewParam.Sort == 1,
        SP = zSortPairs(File,SP,ViewParam);
        ViewParam.Sort = 0;
      else
        for i=1:length(ViewParam.Mode),
          switch ViewParam.Mode(i),
            case 1, FigsDone = zScatterPairs(File,SP,Param,ViewParam);
                    ViewParam.FigNum = ViewParam.FigNum + FigsDone;
            case 2, zListPairs(File,SP,2,ViewParam);
            case 3, zListPairs(File,SP,3,ViewParam);
            case 4, [File,SP,ViewParam]= zDisplayPairs(File,SP,ViewParam);
                    for f=1:length(File),
                      if File(f).Modified == 1,
                        zWriteHandFile(File(f));
                        File(f).Modified = 0;
                        zSaveNTData(File(f));
                      end
                    end
            case 5, zContextViewer(File,SP,Param,ViewParam);
            case 6, 
  CL = zClassLimits;
  for k=1:length(SP),
    f  = SP(k).Filenum;
    p  = File(f).Pair(SP(k).PairIndex);        % Current pair
    N1 = File(f).NT(p.Base2Index);             % reverse order of nucleotides
    N2 = File(f).NT(p.Base1Index);
    sh = (N2.Fit(1,:)-N1.Fit(1,:)) * N1.Rot;   % vector shift from 1 to 2
    p2 = zAnalyzePairFast(N1,N2,CL,sh);        % analyze and classify pair
    F.Pair(k) = p2;                            % make a new file with pairs
    sp2(k) = SP(k);
    sp2(k).Filenum = 1;
    sp2(k).PairIndex = k;
  end
  FigsDone = zScatterPairs(F,sp2,Param,ViewParam);
  ViewParam.FigNum = ViewParam.FigNum + FigsDone;
          end
        end
        if nargin == 3,
          ViewParam.Mode = 0;
        end
      end
    end
  end

  if nargin == 3,
    ViewParam.FigNum = 0;
  end
end  
