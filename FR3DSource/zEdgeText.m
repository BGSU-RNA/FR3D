% zEdgeText converts internal codes for pair interactions into text

% e is the internal code
% Detail = 1 says to give information about subcategories
% C1 and C2 are the codes (1-4) of the two bases making the interaction,
% so we know whether to keep the detail of lowercase letters, eg. cWw
% If C1 is present and C2 is missing, C1 is interpreted as a paircode

% When length(e) == 1, it returns a character string
% When length(e) ~= 1, it returns an empty character string and displays a note

function [E] = zEdgeText(e,Detail,C1,C2)

if nargin < 1,
  e = 1;
end

if nargin < 2,
  Detail = 0;
end

E = [];
EE = [];

if length(e) == 1,
  i = 1;
  switch fix(e(i))
    case    1,    E = [E 'cWw']; EE = [EE 'cWW'];
    case    2,    E = [E 'tWw']; EE = [EE 'tWW'];
    case    3,    E = [E 'cWH']; EE = [EE 'cWH'];
    case    4,    E = [E 'tWH']; EE = [EE 'tWH'];
    case    5,    E = [E 'cWS']; EE = [EE 'cWS'];
    case    6,    E = [E 'tWS']; EE = [EE 'tWS'];
    case    7,    E = [E 'cHh']; EE = [EE 'cHH'];
    case    8,    E = [E 'tHh']; EE = [EE 'tHH'];
    case    9,    E = [E 'cHS']; EE = [EE 'cHS'];
    case   10,    E = [E 'tHS']; EE = [EE 'tHS'];
    case   11,    E = [E 'csS']; EE = [EE 'cSS'];
    case   12,    E = [E 'tSs']; EE = [EE 'tSS'];
    case   13,    E = [E 'bif']; EE = [EE 'bif'];
    case   14,    E = [E 'wat']; EE = [EE 'wat'];
    case   15,    E = [E 'rib']; EE = [EE 'rib'];
    case   -1,    E = [E 'cwW']; EE = [EE 'cWW'];
    case   -2,    E = [E 'twW']; EE = [EE 'tWW'];
    case   -3,    E = [E 'cHW']; EE = [EE 'cHW'];
    case   -4,    E = [E 'tHW']; EE = [EE 'tHW'];
    case   -5,    E = [E 'cSW']; EE = [EE 'cSW'];
    case   -6,    E = [E 'tSW']; EE = [EE 'tSW'];
    case   -7,    E = [E 'chH']; EE = [EE 'cHH'];
    case   -8,    E = [E 'thH']; EE = [EE 'tHH'];
    case   -9,    E = [E 'cSH']; EE = [EE 'cSH'];
    case  -10,    E = [E 'tSH']; EE = [EE 'tSH'];
    case  -11,    E = [E 'cSs']; EE = [EE 'cSS'];
    case  -12,    E = [E 'tsS']; EE = [EE 'tSS'];
    case  -13,    E = [E 'bif']; EE = [EE 'bif'];
    case  -14,    E = [E 'wat']; EE = [EE 'wat'];
    case  -15,    E = [E 'rib']; EE = [EE 'rib'];
    case   21,    E = [E 's35']; EE = [EE 's35'];
    case  -21,    E = [E 's53']; EE = [EE 's53'];
    case   22,    E = [E 's33']; EE = [EE 's33'];
    case  -22,    E = [E 's33']; EE = [EE 's33'];
    case   23,    E = [E 's55']; EE = [EE 's55'];
    case  -23,    E = [E 's55']; EE = [EE 's55'];
    case   28,    E = [E 'perp']; EE = [EE 'perp'];
    case  -28,    E = [E 'perp']; EE = [EE 'perp'];
    case  101,    E = [E 'ncWw']; EE = [EE 'ncWW'];
    case  102,    E = [E 'ntWw']; EE = [EE 'ntWW'];
    case  103,    E = [E 'ncWH']; EE = [EE 'ncWH'];
    case  104,    E = [E 'ntWH']; EE = [EE 'ntWH'];
    case  105,    E = [E 'ncWS']; EE = [EE 'ncWS'];
    case  106,    E = [E 'ntWS']; EE = [EE 'ntWS'];
    case  107,    E = [E 'ncHh']; EE = [EE 'ncHH'];
    case  108,    E = [E 'ntHh']; EE = [EE 'ntHH'];
    case  109,    E = [E 'ncHS']; EE = [EE 'ncHS'];
    case  110,    E = [E 'ntHS']; EE = [EE 'ntHS'];
    case  111,    E = [E 'ncsS']; EE = [EE 'ncSS'];
    case  112,    E = [E 'ntSs']; EE = [EE 'ntSS'];
    case  113,    E = [E 'nbif']; EE = [EE 'nbif'];
    case  114,    E = [E 'nRib']; EE = [EE 'nRib'];
    case -101,    E = [E 'ncwW']; EE = [EE 'ncWW'];
    case -102,    E = [E 'ntwW']; EE = [EE 'ntWW'];
    case -103,    E = [E 'ncHW']; EE = [EE 'ncHW'];
    case -104,    E = [E 'ntHW']; EE = [EE 'ntHW'];
    case -105,    E = [E 'ncSW']; EE = [EE 'ncSW'];
    case -106,    E = [E 'ntSW']; EE = [EE 'ntSW'];
    case -107,    E = [E 'nchH']; EE = [EE 'ncHH'];
    case -108,    E = [E 'nthH']; EE = [EE 'ntHH'];
    case -109,    E = [E 'ncSH']; EE = [EE 'ncSH'];
    case -110,    E = [E 'ntSH']; EE = [EE 'ntSH'];
    case -111,    E = [E 'ncSs']; EE = [EE 'ncSS'];
    case -112,    E = [E 'ntsS']; EE = [EE 'ntSS'];
    case -113,    E = [E 'nbif']; EE = [EE 'nbif'];
    case -114,    E = [E 'nRib']; EE = [EE 'nRib'];
    case  121,    E = [E 'ns35']; EE = [EE 'ns35'];
    case -121,    E = [E 'ns53']; EE = [EE 'ns53'];
    case  122,    E = [E 'ns33']; EE = [EE 'ns33'];
    case -122,    E = [E 'ns33']; EE = [EE 'ns33'];
    case  123,    E = [E 'ns55']; EE = [EE 'ns55'];
    case -123,    E = [E 'ns55']; EE = [EE 'ns55'];
    otherwise     E = [E '----']; EE = [EE '----'];
  end

  if nargin == 4,
    if any(fix(abs(e)) == [11 12 111 112]),
      EE = E;                                     % reset to original
    elseif any(abs(e) == [1 7 8]) && (C1 == C2),
      EE = E;
    end
  end

  if nargin == 3,
    if any(fix(abs(e)) == [11 12 111 112]),
      EE = E;                                     % reset to original
    elseif any(abs(e) == [1 7 8]) && (any(C1 == [1 6 11 16])),
      EE = E;
    end
  end

  if any(fix(abs(e)) == [13 14 28 113 114]),
    EE = E;
  end

  E = EE;

  % add subcategories to basepairs, like cWWa

  if abs(e(i)) < 29 || Detail > 1,                     % pairing or stacking
    if (fix(e(i)) == e(i)) || (Detail == 0),
      E = [E ' '];
    else
      d = abs(e(i)) - fix(abs(e(i)));
      d = fix(d*10+0.0001);
      d = max(1,d);
      T = 'abcdefghijkl';
      E = [E T(d)];
    end
  else
    E = [E ' '];                                       % near pairs will align better now
  end
else
  E = '';
  disp(['zEdgeText requires a single interaction code to be passed in; what was passed in had length ' num2str(length(e))]);
end

