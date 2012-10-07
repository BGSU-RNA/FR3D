% zIsoDiscrepancy(M1,M2,N1,N2) calculates the isodiscrepancy between the pairs represented by nucleotides M1,M2 and N1,N2

function [d,ang,t1,t2,cp,angdiff1,angdiff2,od] = zIsoDiscrepancy(M1,M2,N1,N2)

    R1 = M2.Rot' * M1.Rot;
    R2 = N2.Rot' * N1.Rot;

    ang = zAngleOfRotation(R1 * R2');

    t1 = (M2.Sugar(1,:) - M1.Sugar(1,:))*M1.Rot - (N2.Sugar(1,:) - N1.Sugar(1,:))*N1.Rot;
    t2 = (M1.Sugar(1,:) - M2.Sugar(1,:))*M2.Rot - (N1.Sugar(1,:) - N2.Sugar(1,:))*N2.Rot;

    cp = abs(norm(M2.Sugar(1,:) - M1.Sugar(1,:)) ...
           - norm(N2.Sugar(1,:) - N1.Sugar(1,:)));

    gly1 = (M2.Fit(1,:) - M2.Sugar(1,:))*M1.Rot;
    a1   = atan2(gly1(2),gly1(1));

    gly2 = (N2.Fit(1,:) - N2.Sugar(1,:))*N1.Rot;
    a2   = atan2(gly2(2),gly2(1));

    gly3 = (M1.Fit(1,:) - M1.Sugar(1,:))*M2.Rot;
    a3   = atan2(gly3(2),gly3(1));

    gly4 = (N1.Fit(1,:) - N1.Sugar(1,:))*N2.Rot;
    a4   = atan2(gly4(2),gly4(1));

    angdiff1 = a1 - a2;
    angdiff2 = a3 - a4;

    angdiff1 = angdiff1 - 2*pi*round(angdiff1/(2*pi)); % subtract nearest mult of 360
    angdiff2 = angdiff2 - 2*pi*round(angdiff2/(2*pi)); % subtract nearest mult of 360

    % calculate isodiscrepancy

    od = sqrt((4*ang)^2 + (t1*t1' + t2*t2')/2 + (3*cp)^2); % original
%    d = sqrt(2*(ang)^2 + 2*(t1*t1' + t2*t2') + 2*(cp)^2);
%    d = sqrt(2*(ang)^2 + 2*(t1(1)^2 + t1(2)^2 + t2(1)^2 + t2(2)^2) + 2*(cp)^2);

    d = cp^2 + (t1(1)^2 + t1(2)^2 + t2(1)^2 + t2(2)^2)/2;

    if sign(R1(3,3)) == sign(R2(3,3)),    % normals point in the same direction
      d = d + 4*(angdiff1^2+angdiff2^2)/2;   % use angle in the plane
    else
      d = d + 9*(pi^2);                      % angle of a flip
    end

    d = sqrt(d);

return

%if (angdiff1 > 150) || (angdiff2 > 150),

[a1 a2 a3 a4 angdiff1 angdiff2]

    File.NT(1) = N1;
    File.NT(2) = N2;
    Param.AtOrigin = 1;
    Param.Sugar = 1;
    clf
    zDisplayNT(File,[1 2],Param);
    view(2)
    pause
%end