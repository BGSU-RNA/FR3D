% xGetEdgeNums parses the entries in the GUI concerning base-base interactions

function [ReqEdge,ExEdge,OKPairs,ExPairs,BP1,BP2,EBP1,EBP2,Flank,Range] = xGetEdgeNums(str)

ReqEdge = [];
ExEdge  = [];
OKPairs  = [];
ExPairs  = [];
BP1      = [];
BP2      = [];
EBP1     = [];
EBP2     = [];
Flank    = [];
Range    = [];
RRange   = [];

% ------------------------------ Define relevant strings and associated codes

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

EdgeStr{35} = 'cWw';
BPequiv{35} = [1];

EdgeStr{36} = 'cwW';
BPequiv{36} = [-1];

EdgeStr{37} = 'cHh';
BPequiv{37} = [7];

EdgeStr{38} = 'chH';
BPequiv{38} = [-7];

EdgeStr{39} = 'tHh';
BPequiv{39} = [8];

EdgeStr{40} = 'thH';
BPequiv{40} = [-8];


BPStr{1}    = 'BP';
basephoscode{1}  = 1:17;

BPStr{2}    = '1BP';
basephoscode{2}  = [10];

BPStr{3}    = '2BP';
basephoscode{3}  = [1];

BPStr{4}    = '3BP';
basephoscode{4}  = [11];

BPStr{5}    = '4BP';
basephoscode{5}  = [12 15];

BPStr{6}    = '5BP';
basephoscode{6}  = [2 5 13];

BPStr{7}    = '6BP';
basephoscode{7}  = [3 6];

BPStr{8}    = '7BP';
basephoscode{8}  = [7];

BPStr{9}    = '8BP';
basephoscode{9}  = [8];

BPStr{10}    = '9BP';
basephoscode{10}  = [4 9 14 17];

BPStr{21}    = 'PB';
basephoscode{21}  = 1:17;

BPStr{22}    = '1PB';
basephoscode{22}  = [10];

BPStr{23}    = '2PB';
basephoscode{23}  = [1];

BPStr{24}    = '3PB';
basephoscode{24}  = [11];

BPStr{25}    = '4PB';
basephoscode{25}  = [12 15];

BPStr{26}    = '5PB';
basephoscode{26}  = [2 5 13];

BPStr{27}    = '6PB';
basephoscode{27}  = [3 6];

BPStr{28}    = '7PB';
basephoscode{28}  = [7];

BPStr{29}    = '8PB';
basephoscode{29}  = [8];

BPStr{30}    = '9PB';
basephoscode{30}  = [4 9 14 17];

BPequiv{100} = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

Pairs = {'AA' 'CA' 'GA' 'UA' 'AC' 'CC' 'GC' 'UC' 'AG' 'CG' 'GG' 'UG' 'AU' 'CU' 'GU' 'UU' ''};

str    = strrep(str,'any','');      % no restriction
  
str    = strrep(str,';',',');       % replace delims by commas
str    = strrep(str,':',',');       % replace delims by commas
str    = strrep(str,'|',',');       % replace delims by commas
str    = strrep(str,' ',',');       % replace delims by commas

while strfind(str,',,'),
  str = strrep(str,',,',',');
end

if str(end) == ',',
  str = str(1:(end-1));
end

if str(1) == ',',
  str = str(2:end);
end

commas = strfind(str,',');               % find locations of commas

if isempty(str),
  lim = [];
else
  lim = [0 commas length(str)+1];        % locations of edge specifications
end

for i=1:length(lim)-1                    % loop through tokens
  Token = str(lim(i)+1:lim(i+1)-1);     % extract next token
  if Token(1) == '~',                   % this token, opposite of usual sense
    Token = Token(2:length(Token));
    Reverse = 1;
  else
    Reverse = 0;
  end

  if ((Token(1) == 'n') | (Token(1) == 'N')) && (length(Token) > 1),  % near
    Token = Token(2:length(Token));
    Near = 1;
  else
    Near = 0;
  end

  EdgeNum  = str2num(Token);               % try converting to a number
  PairCode = [];                           % default; nothing
  newBP1   = [];
  newBP2   = [];

  if isempty(EdgeNum)                      % Token IS a string
    edg  = find(strcmp(EdgeStr,Token));    % case sensitive
    edgi = find(strcmpi(EdgeStr,Token));   % case insensitive
    pair = find(strcmpi(Pairs,Token));     % case insensitive
    basephos = find(strcmpi(BPStr,Token)); % case insensitive

    EdgeCode = 99;                         % default; nothing

    if ~isempty(edg),
      EdgeCode = edg(1);
    elseif ~isempty(edgi),
      EdgeCode = edgi(1);
    elseif ~isempty(pair),
      PairCode = pair(1);
    elseif ~isempty(basephos),
      if basephos(1) < 20,
        newBP1 = basephoscode{basephos(1)};
      else
        newBP2 = basephoscode{basephos(1)};          
      end
    end

    newRange = [];

    if strcmpi(Token,'flank') || strcmpi(Token,'f'),
      Flank = 1 - Reverse;
    elseif strcmpi(Token,'local') || strcmpi(Token,'L'),
      newRange = [1 10];
      newRRange = [11 Inf];
    elseif strcmpi(Token,'ested') || strcmpi(Token,'n'),
      newRange = [0 0];
      newRRange = [1 Inf];
    elseif strcmpi(Token,'long-range') || strcmpi(Token,'LR') || strcmpi(Token,'D') || strcmpi(Token,'Distant'),
      newRange = [11 Inf];
      newRRange = [0 10];
    end

    if isempty(Range),
      if ~isempty(newRange),
        Range = newRange;
        RRange = newRRange;
      end
    elseif ~isempty(newRange),
      Range(1) = min(Range(1),newRange(1));     % expand range
      Range(2) = max(Range(2),newRange(2));
      RRange(1) = max(RRange(1),newRRange(1));  % narrow RRange
      RRange(2) = min(RRange(2),newRRange(2));
    end

    EdgeNum = BPequiv{EdgeCode};

    if Near > 0,
      NearEdgeNum = [];
      for j = 1:length(EdgeNum),
        e = EdgeNum(j);
        NearEdgeNum = [NearEdgeNum sign(e) * (100 + abs(e))];
      end
      EdgeNum = NearEdgeNum;

      newBP1 = newBP1 + 100;
      newBP2 = newBP2 + 100;
    end
  end

  if Reverse == 0,
    ReqEdge = [ReqEdge EdgeNum];
    OKPairs = [OKPairs PairCode];
    BP1     = [BP1 newBP1];
    BP2     = [BP2 newBP2];
  else
    ExEdge  = [ExEdge EdgeNum];
    ExPairs = [ExPairs PairCode];
    EBP1    = [EBP1 newBP1];
    EBP2    = [EBP2 newBP2];
    Range   = RRange;
  end
end
