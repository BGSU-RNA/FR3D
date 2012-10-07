% zClassifyPair(N1,N2) calculates the rotation matrix, axis, angle, and shift
% between bases in File that are close enough to possibly be interacting, then
% classifies the interaction

function [Pair,s,coplanar] = zClassifyPair(N1,N2,CL,Exemplar,Force,Verbose)

if nargin < 5,
  Force = 0;
end

if nargin < 6,
  Verbose = 1;
end

if nargin < 3,
  CL = zClassLimits;                              % read ClassLimits matrix
end

if nargin < 4,
  if exist('PairExemplars.mat','file') > 0,
    load('PairExemplars','Exemplar');
  else
    Exemplar = [];
  end
end

% Paircode list
% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

  paircode = 4*(N2.Code-1) + N1.Code;           % AA is 1, CA is 2, etc.

  if Force == 0,
   switch paircode
    case {2, 3, 4, 8, 10, 12},                  % put N2 at the origin
      M1 = N2;
      M2 = N1;
      s  = -1;                                  % bases in reversed order
    otherwise
      M1 = N1;
      M2 = N2;
      s  = 1;                                   % bases in original order
   end
  else
   M1 = N1;
   M2 = N2;
   s  = 1;
  end

  sh = (M2.Fit(1,:)-M1.Fit(1,:)) * M1.Rot;   % vector shift from 1 to 2
                                             % between glycosidic atoms,
                                             % relative to the plane of base 1

  coplanar = 0;                              % default value

  if (abs(sh(3)) < 5) || (Force > 0)
                                             % if small vertical shift
    Pair = zAnalyzePair(M1,M2,CL,Exemplar,sh,Verbose); % analyze and classify pair

    coplanar = Pair.Coplanar;

    if (abs(Pair.Class) >= 30) && (M1.Code == M2.Code) && (Force == 0),  % re-analyze AA CC ...
      M2 = N1;                               % reverse roles of nucleotides
      M1 = N2;
      s  = -1;
      sh2 = (M2.Fit(1,:)-M1.Fit(1,:)) * M1.Rot;   % vector shift from 1 to 2
      Pair2 = zAnalyzePair(M1,M2,CL,Exemplar,sh2,Verbose);%put other base at origin
      if fix(abs(Pair2.Class)) ~= 30,        % some known interaction, or near
        Pair = Pair2;                        % matched with M2 at origin
      else                                   % other interaction, like stacking
        Pair.Class = Pair2.Class;            % original order, use 2nd class
        s = 1;
      end
    end

    if (fix(abs(Pair.Edge)) == 30) && (Force == 0),% remove unclassified pairs
      Pair = [];
    end

  else
    Pair = [];
  end   % if small vertical distance or hand classified
