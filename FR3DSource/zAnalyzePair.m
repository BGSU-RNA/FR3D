% zAnalyzePair(N1,N2,CL) computes distances, angles, and classification
% codes.

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

  a = zCheckCutoffs(Pair.Displ,Pair.Normal,Pair.Ang,Pair.Gap,CL(:,:,Pair.Paircode));
                                           % find possible classifications

  % ------------------------ check for coplanarity

  Pair.Coplanar = 0;                      % default value, not coplanar

  % Criteria for being coplanar or near coplanar:
  %   Pair.Gap must be < 97th percentile among basepairs (1.5179 Angstroms)
  %   Pair.MinDist must be < 97th percentile among basepairs (2.4589 A)
  %   Angle between center-center vector and normals must be > 70.2388 degrees
  %   Angle between normal vectors must be < 39.1315 degrees

  if (abs(Pair.Gap) < 1.5179) && (Pair.MinDist < 2.4589),
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
        dot3Val = 1+(-dot3-(-0.9509))*(-7.4217); % Anton 7/15/2011. For compatibility with octave
%         dot3Val = 1+(-dot3--0.9509)*(-7.4217);        
      elseif -dot3 < -0.7757,           % 97th percentile
        dot3Val = 0.5+(-dot3-(-0.8835))*(-4.6390);  % Anton 7/15/2011. For compatibility with octave
%         dot3Val = 0.5+(-dot3--0.8835)*(-4.6390);
      else
        dot3Val = 0;
      end

      if Pair.MinDist <  1.8982,               % 70th percentile
        MinDistVal = 1;
      elseif Pair.MinDist <  2.1357,           % 90th percentile
        MinDistVal = 1+(Pair.MinDist- 1.8982)*(-2.1050);
      elseif Pair.MinDist <  2.4859,           % 97th percentile
        MinDistVal = 0.5+(Pair.MinDist- 2.1357)*(-1.4280);
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

  % ---------- Notify and remove multiple classifications

  if length(a) > 1,
    if max(fix(a)) > min(fix(a)),              % different integer parts
      if Verbose > 1,
        fprintf('Bases %1s%5s(%1s) and %1s%5s(%1s) fall into categories ', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain);
        for k=1:length(a),
          fprintf('%6.2f ',a(k));
        end
        fprintf('\n');
      end
      a = a(1);
    else
      a = sign(a(1))*min(abs(a));               % use primary version of class
    end
  end

  if Verbose > 2,
    fprintf('Classification is %7.2f after checking cutoffs.\n', a);
  end

  % ---------- Calculate hydrogen bonds for base pairing interactions

  if (abs(a)<14) && (abs(a) - fix(abs(a)) < 0.5),% standard planar interaction
    Pair.Hydrogen = zCheckHydrogen(N1,N2,fix(a));
  else
    Pair.Hydrogen = [];
  end

  % ---------- Eliminate out of plane interactions using Gap cutoffs

  if (abs(a) < 11) || ((abs(a) >= 13) && (abs(a) < 14)),   % avoid cSS, tSS
    if length(Pair.Hydrogen) > 0,
      if abs(Pair.Gap) > 2.0,
        a = 30.1;
      end
    elseif abs(Pair.Gap) > 1.6,
      a = 30.2;
    end
  end

  if Verbose > 2,
    fprintf('Classification is %7.2f after gap cutoffs.\n', a);
  end

  % ---------- Eliminate bad hydrogen bonds

  if (abs(a)<14) && (abs(a) - fix(abs(a)) < 0.5),% standard planar interaction
    if (length(Pair.Hydrogen) > 0),             % there are h-bonds to check
      goodhydrogens = 0;                        % count good bonds
      for h = 1:length(Pair.Hydrogen),
        if isempty(Pair.Hydrogen(h).Angle),     % no hydrogen present
          if Pair.Hydrogen(h).Distance <= 4.5,  % heavy-heavy length < 4.5
            goodhydrogens = goodhydrogens + 1;
          end
        else                                    % hydrogen is present
          if (Pair.Hydrogen(h).Angle >= 110)&&(Pair.Hydrogen(h).Distance <= 4),
            goodhydrogens = goodhydrogens + 1;
          end
        end
      end

      if goodhydrogens < length(Pair.Hydrogen),  % missing hydrogen bond
        gh = length(Pair.Hydrogen);              % required number

        if (length(Pair.Hydrogen) == 4),         % 3 out of 4 is good enough
          gh = 3;
        elseif (length(Pair.Hydrogen) == 3),     % 2 out of 3 is good enough
          gh = 2;
        end

        % Implement a few exceptions to hydrogen bonding rules

        if (fix(abs(a)) == 13) && (Pair.Paircode == 1)
          gh = 0;
        elseif (fix(a) == 13) && (Pair.Paircode == 5)
          gh = 0;
        elseif (fix(a) == -13) && (Pair.Paircode == 5)
          gh = 1;
        elseif (fix(abs(a)) == 13) && (Pair.Paircode == 6)
          gh = 1;
        elseif (fix(a) == 11) && (Pair.Paircode == 6)
          gh = 0;
        elseif (fix(a) == 12) && (Pair.Paircode == 11)      % GG tSS
          gh = 2;
        elseif (a == 9) && (Pair.Paircode == 6)              % CC cHS
          gh = 1;
        end
     
        if goodhydrogens < gh,
          a = 30.3;                              % reject this pair
        end
      end
    end
  end

  if Verbose > 2,
    fprintf('Classification is %7.2f.  %d good hydrogens required, %d found.\n', a, gh, goodhydrogens);
  end

  % ------------ If not base pairing, check for stacking

  Pair.StackingOverlap = 0;

  if fix(a) == 30,

    % ---------- Measure stacking overlap

    SO1 = zStackingOverlap(N1,N2);
    SO2 = zStackingOverlap(N2,N1);

    if (SO1 > 0) && (SO2 > 0),
      Pair.StackingOverlap = (SO1+SO2)/2;
    else
      Pair.StackingOverlap = -max(SO1,SO2);
    end

    % ---------- Provisional stacking category  

    if Pair.Displ(3) > 0,
      if Pair.Normal(3) > 0,
        aa = 21;                       % second base above, pointing up
        aaa = 121;
      else
        aa = 22;                       % second base above, pointing down
        aaa = 122;
      end
    else
      if Pair.Normal(3) > 0,
        aa = -21;                      % second base below, pointing up
        aaa = -121;
      else
        aa = 23;                       % second base below, pointing down
        aaa = 123;
      end
    end

    % ----------- Is an unclassified pair really stacked?

    cond(1) = Pair.StackingOverlap > 0;
    cond(2) = Pair.MinDist < 4;
    cond(3) = Pair.MinDist > 1;
    cond(4) = abs(Pair.Normal(3)) > 0.6;
    cond(5) = Pair.StackingOverlap < 0;
    cond(6) = Pair.MinDist < 4.5;
    cond(7) = abs(Pair.Normal(3)) > 0.5;

    if sum(cond(1:4)) == 4,
       a = aa;
    elseif (cond(1)+cond(5) > 0) && (cond(6) > 0) && (cond(3) > 0) && cond(7) > 0,
       a = aaa;
    end          
  end

  % ----------- If classified as stacking, is there really overlap?

%  if (fix(a) >= 21) && (fix(a) < 24) && (Pair.StackingOverlap == 0),
%    a = 30;
%    if N1.Loc(1,1) < N2.Loc(1,1),     % only print each pair once
      % fprintf('Bases %1s%5s(%1s) and %1s%5s(%1s) have no stacking overlap\n', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain);
%    end
%  end

  % ----------------------- Perpendicular interactions

  if fix(a) == 30,                             % not yet classified
    if (Pair.Displ(1) > -6.0) && ...
       (Pair.Displ(1) <  6.0) && ...
       (Pair.Displ(2) > -6.0) && ...
       (Pair.Displ(2) <  6.0) && ...
       (Pair.Displ(3) > -9.0) && ...
       (Pair.Displ(3) <  9.0) && ...
       (abs(Pair.Normal(3)) <= 0.5) && ...
       (Pair.MinDist < 4),
      a = 28;
    end
  end

  % ----------------- Record nearest basepair exemplar if no classification yet

  if (fix(a) == 30) && abs(Pair.Gap) < 2.25 && Pair.MinDist < 4.75,

    % -------------------------- find distance to nearest exemplar

    if ~isempty(Exemplar),
      [c,d,ff,gg,h] = zDistanceToExemplars(Exemplar,N1,N2);
      Pair.Classes   = c(1:3);
      Pair.Distances = d(1:3);
    else
      Pair.Classes   = 99 * ones(1,3);
      Pair.Distances = 99999999 * ones(1,3);
    end

    % --------------------------- check hydrogen bonds if nearest exemplar
    % --------- is a basepair, even if this pair is far from that exemplar

    if (abs(a) >= 14) && (abs(Pair.Classes(1)) < 14),
      aa = abs(Pair.Classes(1));
      if aa - fix(aa) < 0.5,                       % usual hydrogen bonds
        Pair.Hydrogen = zCheckHydrogen(N1,N2,fix(Pair.Classes(1)));
      end
    end

    if (Pair.Distances(1) < 0.8) && (Pair.Normal(3) * Exemplar(ff(1),gg(1)).R(3,3) > 0), % same flip
      b = a-fix(a);                          % extract decimal code for reason
      c = Pair.Classes(1);
      a = sign(c) * (100 + abs(c) + b/1000);
    elseif (Pair.Distances(2) < 0.8) && (Pair.Normal(3) * Exemplar(ff(2),gg(2)).R(3,3) > 0), % same flip
      b = a-fix(a);                          % extract decimal code for reason
      c = Pair.Classes(2);
      a = sign(c) * (100 + abs(c) + b/1000);
    elseif (Pair.Distances(3) < 0.8) && (Pair.Normal(3) * Exemplar(ff(3),gg(3)).R(3,3) > 0), % same flip
      b = a-fix(a);                          % extract decimal code for reason
      c = Pair.Classes(3);
      a = sign(c) * (100 + abs(c) + b/1000);
    end
  else
    Pair.Classes   = 99 * ones(1,3);
    Pair.Distances = 99999999 * ones(1,3);
  end

  % ------------------------ store the classification

  Pair.Class = a;

  % ------------------------ store the edge information

  % reverse classification for GC and CG pairs, but why, exactly????

  if ((Pair.Paircode == 7) || (Pair.Paircode == 10)) && (mod(abs(Pair.Class),100) < 14),
    Pair.Edge = -Pair.Class;
  else
    Pair.Edge = Pair.Class;
  end

  Pair.EdgeText = zEdgeText(Pair.Edge,1,Pair.Paircode);  % ever used?


if 0 > 1,
  File = zAddNTData('3i8i');
  load('PairExemplars','Exemplar');
  CL = zClassLimits;  
  N1 = File.NT(zIndexLookup(File,'1981'));
  N2 = File.NT(zIndexLookup(File,'758'));
  Displ = (N2.Fit(1,:)-N1.Fit(1,:)) * N1.Rot;   % vector shift from 1 to 2
  p = zAnalyzePair(N1,N2,CL,Exemplar,Displ,3)
end
