% xGetEdgeNums parses the entries in the GUI concerning base-base interactions

function [ReqEdge,ExEdge,OKPairs,ExPairs,BP1,BP2,EBP1,EBP2,BR1,BR2,EBR1,EBR2,Flank,Range,Coplanar,ReqBB,ExBB] = xGetEdgeNums(str)

ReqEdge = [];
ExEdge  = [];
OKPairs  = [];
ExPairs  = [];
BP1      = [];
BP2      = [];
EBP1     = [];
EBP2     = [];
BR1      = [];
BR2      = [];
EBR1     = [];
EBR2     = [];
Flank    = [];
Range    = [];
Coplanar = 0;
NearCoplanar = 0;
ReqBB    = [];
ExBB     = [];

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

EdgeStr{41} = 'perp';
BPequiv{41} = [28 -28];

EdgeStr{42} = 'perpendicular';
BPequiv{42} = [28 -28];

EdgeStr{43} = 'wat';
BPequiv{43} = [14 -14];

%BPCat = [2 6 8 0 6 7 8 9 0 1 3 4 5 0 5 8 0];  % updated 8-7-2008
%         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 

% Convert external category names for BP interactions to internal codes

BPStr{1}    = 'BPh';
basephoscode{1}  = 1:19;

BPStr{2}    = '1BPh';
basephoscode{2}  = [10];

BPStr{3}    = '2BPh';
basephoscode{3}  = [1];

BPStr{4}    = '3BPh';
basephoscode{4}  = [11];

BPStr{5}    = '4BPh';
basephoscode{5}  = [12 19];

BPStr{6}    = '5BPh';
basephoscode{6}  = [13 15];

BPStr{7}    = '6BPh';
basephoscode{7}  = [2 5];

BPStr{8}    = '7BPh';
basephoscode{8}  = [3 6 18];

BPStr{9}    = '8BPh';
basephoscode{9}  = [7];

BPStr{10}    = '9BPh';
basephoscode{10}  = [8 16];

BPStr{11}    = '0BPh';
basephoscode{11}  = [4 9 14 17];

BPStr{12}    = '4bBP';
basephoscode{12}  = [19];

BPStr{13}    = '8bBP';
basephoscode{13}  = [18];

BPStr{14}    = '4bBPh';
basephoscode{14}  = [19];

BPStr{15}    = '8bBPh';
basephoscode{15}  = [18];

BPStr{21}    = 'PhB';
basephoscode{21}  = 1:19;

BPStr{22}    = '1PhB';
basephoscode{22}  = [10];

BPStr{23}    = '2PhB';
basephoscode{23}  = [1];

BPStr{24}    = '3PhB';
basephoscode{24}  = [11];

BPStr{25}    = '4PhB';
basephoscode{25}  = [12];

BPStr{26}    = '5PhB';
basephoscode{26}  = [13 15];

BPStr{27}    = '6PhB';
basephoscode{27}  = [2 5];

BPStr{28}    = '7PhB';
basephoscode{28}  = [3 6 18];

BPStr{29}    = '8PhB';
basephoscode{29}  = [7];

BPStr{30}    = '9PhB';
basephoscode{30}  = [8 16];

BPStr{31}    = '0PhB';
basephoscode{31}  = [4 9 14 17];

BPStr{32}    = '4bPB';
basephoscode{32}  = [19];

BPStr{33}    = '8bPB';
basephoscode{33}  = [18];

BPStr{34}    = '4bPhB';
basephoscode{34}  = [19];

BPStr{35}    = '8bPhB';
basephoscode{35}  = [18];

BPequiv{100} = [];

% Convert external category names for BR interactions to internal codes

BRStr{1}    = 'BR';
baseribosecode{1}  = 1:19;

BRStr{2}    = '1BR';
baseribosecode{2}  = [10];

BRStr{3}    = '2BR';
baseribosecode{3}  = [1];

BRStr{4}    = '3BR';
baseribosecode{4}  = [11];

BRStr{5}    = '4BR';
baseribosecode{5}  = [12 19];

BRStr{6}    = '5BR';
baseribosecode{6}  = [13 15];

BRStr{7}    = '6BR';
baseribosecode{7}  = [2 5];

BRStr{8}    = '7BR';
baseribosecode{8}  = [3 6 18];

BRStr{9}    = '8BR';
baseribosecode{9}  = [7];

BRStr{10}    = '9BR';
baseribosecode{10}  = [8 16];

BRStr{11}    = '0BR';
baseribosecode{11}  = [4 9 14 17];

BRStr{12}    = '4bBP';
baseribosecode{12}  = [19];

BRStr{13}    = '8bBP';
baseribosecode{13}  = [18];

BRStr{14}    = '4bBR';
baseribosecode{14}  = [19];

BRStr{15}    = '8bBR';
baseribosecode{15}  = [18];

BRStr{21}    = 'RB';
baseribosecode{21}  = 1:19;

BRStr{22}    = '1RB';
baseribosecode{22}  = [10];

BRStr{23}    = '2RB';
baseribosecode{23}  = [1];

BRStr{24}    = '3RB';
baseribosecode{24}  = [11];

BRStr{25}    = '4RB';
baseribosecode{25}  = [12];

BRStr{26}    = '5RB';
baseribosecode{26}  = [13 15];

BRStr{27}    = '6RB';
baseribosecode{27}  = [2 5];

BRStr{28}    = '7RB';
baseribosecode{28}  = [3 6 18];

BRStr{29}    = '8RB';
baseribosecode{29}  = [7];

BRStr{30}    = '9RB';
baseribosecode{30}  = [8 16];

BRStr{31}    = '0RB';
baseribosecode{31}  = [4 9 14 17];

BRStr{32}    = '4bRB';
baseribosecode{32}  = [19];

BRStr{33}    = '8bRB';
baseribosecode{33}  = [18];

BRStr{34}    = '4bRB';
baseribosecode{34}  = [19];

BRStr{35}    = '8bRB';
baseribosecode{35}  = [18];

BPequiv{100} = [];

zBackboneCodes;
BBStr = Codes;

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
  newBR1   = [];
  newBR2   = [];
  CP       = 0;
  nCP      = 0;
  BBCode   = [];

  if ~isempty(EdgeNum)                      % Token is a number

    if any (fix(abs(EdgeNum)) == [1 2 7 8]),
      EdgeNum = [EdgeNum -EdgeNum];
    end

    if Near > 0,
      NearEdgeNum = [];
      for j = 1:length(EdgeNum),
        e = EdgeNum(j);
        NearEdgeNum = [NearEdgeNum sign(e) * (100 + abs(e))];
      end
      EdgeNum = NearEdgeNum;
    end

  else
    edg  = find(strcmp(EdgeStr,Token));    % case sensitive
    edgi = find(strcmpi(EdgeStr,Token));   % case insensitive
    pair = find(strcmpi(Pairs,Token));     % case insensitive
    basephos = find(strcmpi(BPStr,Token)); % case insensitive
    baseribose = find(strcmpi(BRStr,Token)); % case insensitive
    bb   = find(strcmpi(BBStr,Token));

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
    elseif ~isempty(baseribose),
      if baseribose(1) < 20,
        newBR1 = baseribosecode{baseribose(1)};
      else
        newBR2 = baseribosecode{baseribose(1)};          
      end
    elseif ~isempty(bb),
      BBCode = bb(1);
    end

    newRange = [];

    if strcmpi(Token,'flank') || strcmpi(Token,'f') || strcmpi(Token,'flankss') || strcmpi(Token,'borderss'),
      Flank = 1 - Reverse;
    elseif strcmpi(Token,'coplanar') || strcmpi(Token,'cp'),
      CP = 1 - 2*Reverse;
    elseif strcmpi(Token,'local') || strcmpi(Token,'L'),
      newRange = [0 3];
    elseif strcmpi(Token,'ested') || strcmpi(Token,'n'),
      newRange = [0 0];
    elseif strcmpi(Token,'long-range') || strcmpi(Token,'LR') || strcmpi(Token,'D') || strcmpi(Token,'Distant'),
      newRange = [4 Inf];
    elseif ~isempty(strfind(lower(Token),'crossing')),
      j = strfind(Token,'_');
      newRange(1) = str2num(Token((j(1)+1):(j(2)-1)));
      newRange(2) = str2num(Token((j(2)+1):end));
    end

    if isempty(Range),
      if ~isempty(newRange),
        Range = newRange;
      end
    elseif ~isempty(newRange),
      Range(1) = min(Range(1),newRange(1));     % expand range
      Range(2) = max(Range(2),newRange(2));
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
      newBR1 = newBR1 + 100;
      newBR2 = newBR2 + 100;

      if (CP ~= 0),
        nCP = CP;
        CP  = 0;
      end
    end
  end

  if Reverse == 0,
    ReqEdge = [ReqEdge EdgeNum];
    OKPairs = [OKPairs PairCode];
    BP1     = [BP1 newBP1];
    BP2     = [BP2 newBP2];
    BR1     = [BR1 newBR1];
    BR2     = [BR2 newBR2];
    ReqBB   = [ReqBB BBCode];
  else
    ExEdge  = [ExEdge EdgeNum];
    ExPairs = [ExPairs PairCode];
    EBP1    = [EBP1 newBP1];
    EBP2    = [EBP2 newBP2];
    EBR1    = [EBR1 newBR1];
    EBR2    = [EBR2 newBR2];
    ExBB    = [ExBB BBCode];
  end

  if (CP ~= 0),
    Coplanar = CP;
  end
  if (nCP ~= 0)
    NearCoplanar = nCP;
  end
end

if     (Coplanar == 0) && (NearCoplanar == 0),
  Coplanar = [];
elseif (Coplanar == -1) && (NearCoplanar == 0),
  Coplanar = 0;
elseif (Coplanar == 1) && (NearCoplanar == 0),
  Coplanar = 1;
elseif (Coplanar == 1) && (NearCoplanar == -1),
  Coplanar = 1;
elseif (Coplanar == 1) && (NearCoplanar == 1),
  Coplanar = 2;
elseif (Coplanar == 0) && (NearCoplanar == 1),
  Coplanar = 3;
elseif (Coplanar == -1) && (NearCoplanar == 1),
  Coplanar = 3;
elseif (Coplanar == 0) && (NearCoplanar == -1),
  Coplanar = 4;
elseif (Coplanar == -1) && (NearCoplanar == -1),
  Coplanar = 5;
end
