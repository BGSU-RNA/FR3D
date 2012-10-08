% zLooseCoplanar(N1,N2,CL) computes distances, angles, and classification
% codes and gives a looser classification of "coplanar" than the official one,
% but still useful for building SCFGs

function [Pair] = zAnalyzePair(N1,N2,CL,Exemplar,Displ,Verbose)

  if nargin < 3,
    CL = zClassLimits;                              % read ClassLimits matrix

    if exist('PairExemplars.mat','file') > 0,
      load('PairExemplars','Exemplar');
    else
      Exemplar = [];
    end
  end

  if nargin < 6,
    Verbose = 3;
  end

  if nargin < 5,
    Displ = (N2.Fit(1,:)-N1.Fit(1,:)) * N1.Rot;   % vector shift from 1 to 2
  end

  Pair.Paircode = 4*(N2.Code-1) + N1.Code;     % AA is 1, CA is 2, etc.

  Pair.Displ = Displ;                        % vector shift (from 1 to 2)

  ro = N1.Rot'*N2.Rot;                       % rotation matrix from 1 to 2
  Pair.Normal = ro(:,3)';                    % normal to second plane

  if ro(3,3) > 0,                            % depending on orientation of 2,
    [ax,ang] = zAxisAngle(ro);               % rotation angle without a flip
  else
    [ax,ang] = zAxisAngle(ro*diag([-1 1 -1])); % flip base 2 first
  end

  Pair.Rot      = ro;
  Pair.RotAx    = ax';
  Pair.Ang      = ang;

  Pair.PlaneAng = acos(abs(N1.Rot(:,3)'*N2.Rot(:,3)))*57.29577951308232; 
                                             % angle between planes

  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  d = zDistance(N2.Fit(1:Lim(2,N2.Code),:), N1.Center); 
                                           % distances to base 1 center
  [y,m] = min(d);                          % identify the closest atom
  m = m(1);                                % in case of a tie, use the first
  Pair.Gap = N1.Rot(:,3)'*(N2.Fit(m,:)-N1.Center)';% height above plane of 1

  Pair.MinDist = min(min(zDistance(N1.Fit,N2.Fit)));

  % ------------------------ check for coplanarity

  Pair.Coplanar = 0;                      % default value, not coplanar

  % Criteria for being coplanar or near coplanar:
  %   Pair.Gap must be < 97th percentile among basepairs (1.5179 Angstroms)
  %   Pair.MinDist must be < 97th percentile among basepairs (2.4589 A)
  %   Angle between center-center vector and normals must be > 70.2388 degrees
  %   Angle between normal vectors must be < 39.1315 degrees

  if (abs(Pair.Gap) < 1.4) && (Pair.MinDist < 2.8),
    v  = N1.Center - N2.Center;           % vector from center to center
    v  = v / norm(v);                     % normalize

    dot1 = abs(v * N1.Rot(:,3));          % to calculate angle: v and normal
    dot2 = abs(v * N2.Rot(:,3));
    dot3 = abs(N1.Rot(:,3)' * N2.Rot(:,3));

    if (dot1 < 0.3381) && (dot2 < 0.3381) && (dot3 > 0.7757),

      d = zDistance(N1.Fit(1:Lim(2,N1.Code),:), N2.Center); 
                                           % distances to base 2 center
      [y,m] = min(d);                      % identify the closest atom
      m = m(1);                            % in case of a tie, use the first
      Gap2 = N2.Rot(:,3)'*(N1.Fit(m,:)-N2.Center)';% height above plane of 1

      if abs(Pair.Gap) <  0.5062,               % 70th percentile
        Gap1Val = 1;
      elseif abs(Pair.Gap) <  0.9775,           % 90th percentile
        Gap1Val = 1+(abs(Pair.Gap)- 0.5062)*(-1.0609);
      elseif abs(Pair.Gap) <  1.5179,           % 97th percentile
        Gap1Val = 0.5+(abs(Pair.Gap)- 0.9775)*(-0.9252);
      else
        Gap1Val = 0;
      end

      if abs(Gap2) <  0.5062,               % 70th percentile
        Gap2Val = 1;
      elseif abs(Gap2) <  0.9775,           % 90th percentile
        Gap2Val = 1+(abs(Gap2)- 0.5062)*(-1.0609);
      elseif abs(Gap2) <  1.5179,           % 97th percentile
        Gap2Val = 0.5+(abs(Gap2)- 0.9775)*(-0.9252);
      else
        Gap2Val = 0;
      end

      if dot1 <  0.1139,               % 70th percentile
        dot1Val = 1;
      elseif dot1 <  0.2193,           % 90th percentile
        dot1Val = 1+(dot1- 0.1139)*(-4.7408);
      elseif dot1 <  0.3381,           % 97th percentile
        dot1Val = 0.5+(dot1- 0.2193)*(-4.2103);
      else
        dot1Val = 0;
      end

      if dot2 <  0.1139,               % 70th percentile
        dot2Val = 1;
      elseif dot2 <  0.2193,           % 90th percentile
        dot2Val = 1+(dot2- 0.1139)*(-4.7408);
      elseif dot2 <  0.3381,           % 97th percentile
        dot2Val = 0.5+(dot2- 0.2193)*(-4.2103);
      else
        dot2Val = 0;
      end

      if -dot3 < -0.9509,               % 70th percentile
        dot3Val = 1;
      elseif -dot3 < -0.8835,           % 90th percentile
        dot3Val = 1+(-dot3--0.9509)*(-7.4217);
      elseif -dot3 < -0.7757,           % 97th percentile
        dot3Val = 0.5+(-dot3--0.8835)*(-4.6390);
      else
        dot3Val = 0;
      end

      if Pair.MinDist <  1.8982,               % 70th percentile
        MinDistVal = 1;
      elseif Pair.MinDist <  2.1357,           % 90th percentile
        MinDistVal = 1+(Pair.MinDist- 1.8982)*(-2.1050);
      elseif Pair.MinDist <  2.8,              % 97th percentile
        MinDistVal = 0.5+(Pair.MinDist- 2.1357)*(-0.5/(2.8-2.1357));
      else
        MinDistVal = 0;
      end

      % Pair.Coplanar is 1 if all are within the 70th percentile
      % Pair.Coplanar is 0.5 if all are within the 90th percentile
      % Pair.Coplanar is > 0 if all are within the 97th percentile
      % Between these, it decreases linearly

      Pair.Coplanar = min([Gap1Val Gap2Val dot1Val dot2Val dot3Val MinDistVal]);

if 10 < 1,
  Pair.Coplanar
  clf
  F.NT(1) = N1;
  F.NT(2) = N2;
  F.Filename = '--';
  zDisplayNT(F);
  pause
end

    end
  end

