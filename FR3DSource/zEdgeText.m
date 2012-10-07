% zEdgeText converts internal codes for pair interactions into text

% e is the internal code
% Detail = 1 says to give information about subcategories
% C1 and C2 are the codes (1-4) of the two bases making the interaction,
% so we know whether to keep the detail of lowercase letters, eg. cWw
% If C1 is present and C2 is missing, C1 is interpreted as a paircode

function [E] = zEdgeText(e,Detail,C1,C2)

if nargin < 1,
  e = 1;
end

if nargin < 2,
  Detail = 0;
end

E = [];

for i=1:length(e),
  switch fix(e(i))
    case    1,    E = [E 'cWw'];
    case    2,    E = [E 'tWw'];
    case    3,    E = [E 'cWH'];
    case    4,    E = [E 'tWH'];
    case    5,    E = [E 'cWS'];
    case    6,    E = [E 'tWS'];
    case    7,    E = [E 'cHh'];
    case    8,    E = [E 'tHh'];
    case    9,    E = [E 'cHS'];
    case   10,    E = [E 'tHS'];
    case   11,    E = [E 'cSs'];
    case   12,    E = [E 'tSs'];
    case   13,    E = [E 'bif'];
    case   14,    E = [E 'wat'];
    case   15,    E = [E 'rib'];
    case   -1,    E = [E 'cwW'];
    case   -2,    E = [E 'twW'];
    case   -3,    E = [E 'cHW'];
    case   -4,    E = [E 'tHW'];
    case   -5,    E = [E 'cSW'];
    case   -6,    E = [E 'tSW'];
    case   -7,    E = [E 'chH'];
    case   -8,    E = [E 'thH'];
    case   -9,    E = [E 'cSH'];
    case  -10,    E = [E 'tSH'];
    case  -11,    E = [E 'csS'];
    case  -12,    E = [E 'tsS'];
    case  -13,    E = [E 'bif'];
    case  -14,    E = [E 'wat'];
    case  -15,    E = [E 'rib'];
    case   21,    E = [E 's35'];
    case  -21,    E = [E 's53'];
    case   22,    E = [E 's33'];
    case  -22,    E = [E 's33'];
    case   23,    E = [E 's55'];
    case  -23,    E = [E 's55'];
    case   28,    E = [E 'perp'];
    case  -28,    E = [E 'perp'];
    case  101,    E = [E 'ncWw'];
    case  102,    E = [E 'ntWw'];
    case  103,    E = [E 'ncWH'];
    case  104,    E = [E 'ntWH'];
    case  105,    E = [E 'ncWS'];
    case  106,    E = [E 'ntWS'];
    case  107,    E = [E 'ncHh'];
    case  108,    E = [E 'ntHh'];
    case  109,    E = [E 'ncHS'];
    case  110,    E = [E 'ntHS'];
    case  111,    E = [E 'ncSs'];
    case  112,    E = [E 'ntSs'];
    case  113,    E = [E 'nbif'];
    case  114,    E = [E 'nRib'];
    case -101,    E = [E 'ncwW'];
    case -102,    E = [E 'ntwW'];
    case -103,    E = [E 'ncHW'];
    case -104,    E = [E 'ntHW'];
    case -105,    E = [E 'ncSW'];
    case -106,    E = [E 'ntSW'];
    case -107,    E = [E 'nchH'];
    case -108,    E = [E 'nthH'];
    case -109,    E = [E 'ncSH'];
    case -110,    E = [E 'ntSH'];
    case -111,    E = [E 'ncsS'];
    case -112,    E = [E 'ntsS'];
    case -113,    E = [E 'nbif'];
    case -114,    E = [E 'nRib'];
    case  121,    E = [E 'ns35'];
    case -121,    E = [E 'ns53'];
    case  122,    E = [E 'ns33'];
    case -122,    E = [E 'ns33'];
    case  123,    E = [E 'ns55'];
    case -123,    E = [E 'ns55'];
    otherwise     E = [E '----'];
  end


EE = E;
EE((end-1):end) = upper(E((end-1):end));

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
 end
end

