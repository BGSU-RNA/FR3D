%mGetReqInterFromGUI.m

function [ReqEdge,ExEdge,OKPairs,ExPairs] = xGetEdgeNums(str)

if strmatch(str,'any')        % most Matrix boxes in the GUI will be like this
  ReqEdge = [];
  ExEdge  = [];
  OKPairs  = [];
  ExPairs  = [];
else
  ReqEdge = [];
  ExEdge  = [];
  OKPairs  = [];
  ExPairs  = [];

  EdgeStr{1} = 'cWW';
  BPequiv{1} = [1 -1];

  EdgeStr{2} = 'tWW';
  BPequiv{2} = [2 -2];

  EdgeStr{3} = 'cWH';
  BPequiv{3} = [3];

  EdgeStr{4} = 'cHW';
  BPequiv{4} = [-3];

  EdgeStr{5} = 'tWH';
  BPequiv{5} = [4];

  EdgeStr{6} = 'tHW';
  BPequiv{6} = [-4];

  EdgeStr{7} = 'cWS';
  BPequiv{7} = [5];

  EdgeStr{8} = 'cSW';
  BPequiv{8} = [-5];

  EdgeStr{9} = 'tWS';
  BPequiv{9} = [6];

  EdgeStr{10} = 'tSW';
  BPequiv{10} = [-6];

  EdgeStr{11} = 'cHH';
  BPequiv{11} = [7 -7];

  EdgeStr{12} = 'tHH';
  BPequiv{12} = [8 -8];

  EdgeStr{13} = 'cHS';
  BPequiv{13} = [9];

  EdgeStr{14} = 'cSH';
  BPequiv{14} = [-9];

  EdgeStr{15} = 'tHS';
  BPequiv{15} = [10];

  EdgeStr{16} = 'tSH';
  BPequiv{16} = [-10];

  EdgeStr{17} = 'cSS';
  BPequiv{17} = [11 -11];

  EdgeStr{18} = 'cSs';
  BPequiv{18} = [11];

  EdgeStr{19} = 'csS';
  BPequiv{19} = [-11];

  EdgeStr{20} = 'tSS';
  BPequiv{20} = [12 -12];

  EdgeStr{21} = 'tSs';
  BPequiv{21} = [12];

  EdgeStr{22} = 'tsS';
  BPequiv{22} = [-12];

  EdgeStr{23} = 'sP';
  BPequiv{23} = [21 -21];

  EdgeStr{24} = 'sA';
  BPequiv{24} = [22 -22 23 -23];

  EdgeStr{25} = 's35';
  BPequiv{25} = [21];

  EdgeStr{26} = 's53';
  BPequiv{26} = [-21];

  EdgeStr{27} = 's33';
  BPequiv{27} = [22 -22];

  EdgeStr{28} = 's55';
  BPequiv{28} = [23 -23];

  EdgeStr{29} = 'NP';            % close enough to interact,
                                 % but not classified or near a class
  BPequiv{29} = [30 -30];

  EdgeStr{30} = 'Pair';
  BPequiv{30} = [-12:-1 1:12];

  EdgeStr{31} = 'Stack';
  BPequiv{31} = [21:23 -23:-21];

  EdgeStr{32} = 'cis';
  BPequiv{32} = [1 3 5 7 9 11 -1 -3 -5 -7 -9 -11];

  EdgeStr{33} = 'trans';
  BPequiv{33} = [2 4 6 8 10 12 -2 -4 -6 -8 -10 -12];

  EdgeStr{34} = 'bif';
  BPequiv{34} = [13 -13];

  BPequiv{100} = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

  Pairs = {'AA' 'CA' 'GA' 'UA' 'AC' 'CC' 'GC' 'UC' 'AG' 'CG' 'GG' 'UG' 'AU' 'CU' 'GU' 'UU' ''};
    
  str    = regexprep(str,';| ',',');       % replace delims by commas
  while strfind(str,',,'),
    str = regexprep(str,',,',',');
  end
  commas = strfind(str,',');               % find locations of commas
  lim    = [0 commas length(str)+1];       % locations of edge specifications

  for i=1:length(lim)-1                    % loop through tokens
      
      NewStr = str(lim(i)+1:lim(i+1)-1); % extract next token
      if NewStr(1) == '~',                 % opposite of usual sense
        NewStr = NewStr(2:length(NewStr));
        Reverse = 1;
      else
        Reverse = 0;
      end
      if (NewStr(1) == 'n') | (NewStr(1) == 'N'),  % near
        NewStr = NewStr(2:length(NewStr));
        Near = 1;
      else
        Near = 0;
      end
      EdgeNum = str2num(NewStr);               % convert to number
      if isempty(EdgeNum)                      % NewStr IS a string
        edg  = find(strcmp(EdgeStr,NewStr));   % case sensitive
        edgi = find(strcmpi(EdgeStr,NewStr));  % case insensitive
        pair = find(strcmpi(Pairs,NewStr));    % case insensitive

        EdgeCode = 99;                         % default; nothing
        PairCode = [];                         % default; nothing

        if ~isempty(edg),
          EdgeCode = edg(1);
        elseif ~isempty(edgi),
          EdgeCode = edgi(1);
        elseif ~isempty(pair),
          PairCode = pair(1);
        end

        EdgeNum = BPequiv{EdgeCode};

        if Near > 0,
          NearEdgeNum = [];
          J = length(EdgeNum);
          for j = 1:J,
            e = EdgeNum(j);
            NearEdgeNum = [NearEdgeNum sign(e) * (100 + abs(e))];
          end
          EdgeNum = NearEdgeNum;
        end
      end

      if Reverse == 0,
        ReqEdge = [ReqEdge EdgeNum];
        OKPairs = [OKPairs PairCode];
      else
        ExEdge  = [ExEdge EdgeNum];
        ExPairs = [ExPairs PairCode];
      end

      

    end
end



return



% some original code by Ali Mokdad:

    BPcodes= {'cWW','tWW','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHH','tHH','cHS','cSH','tHS','tSH','cSS','cSs','csS','tSS','tSs','tsS','sP','sA'};
    BPequiv1=[1    ,2    ,3    ,-3   ,4    ,-4   ,5    ,-5   ,6    ,-6   ,7    ,8    ,9    ,-9   ,10   ,-10  ,11   ,11   ,-11  ,12   ,12   ,-12  ,15, 16   ];
    BPequiv2=[-1   ,-2   ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,-7   ,-8   ,0    ,0    ,0    ,0    ,-11  ,0    ,0    ,-12  ,0    ,0    ,17, 18   ];
    BPequiv3=[0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0,   0     ];
    BPequiv4=[0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0,   0     ];

for k=1:length(BPcodes),
  fprintf('  EdgeStr{%d} = ''%s'';\n', k, BPcodes{k});
  fprintf('  BPequiv{%d} = [', k);
  if BPequiv1(k) ~= 0,
    fprintf('%d', BPequiv1(k));
  end
  if BPequiv2(k) ~= 0,
    fprintf(' %d', BPequiv2(k));
  end
  if BPequiv3(k) ~= 0,
    fprintf(' %d', BPequiv3(k));
  end
  if BPequiv4(k) ~= 0,
    fprintf(' %d', BPequiv4(k));
  end
  fprintf('];\n\n');
end
