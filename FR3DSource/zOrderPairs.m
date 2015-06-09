% zOrderPairs calculates discrepancy between the first LMax pairs in List.  It gives a penalty for non-coplanar pairs.  It returns the "best" pairs, but a better method in zFindExemplars supersedes this choice.

function [PD,ID,i,k] = zOrderPairs(File,List,LMax,Verbose)

    L = min(LMax,length(List(:,1))); % Limit the number of pairs to consider

    PD = zeros(L,L);               % Initialize the pair discrepancy
    ID = zeros(L,L);
    Ang = zeros(L,L);
    T1  = zeros(L,L);
    T2  = zeros(L,L);
    CP  = zeros(L,L);
    for k = 1:L,                   % Slow nested loop
      f1    = List(k,3);
      Model = List(k,[1 2]);
      NT1 = File(f1).NT(Model(1));
      NT2 = File(f1).NT(Model(2));

      if File(f1).Coplanar(Model(1),Model(2)) == 0,
        PD(k,k) = 900;
        ID(k,k) = 900;
      elseif File(f1).Coplanar(Model(1),Model(2)) < 0.5,
%        PD(k,k) = 100;
%        ID(k,k) = 100;
      end

      for m = (k+1):L,
        f2    = List(m,3);
        Cand  = List(m,[1 2]);
        PD(k,m) = xDiscrepancy(File(f1),Model,File(f2),Cand);

        NT3 = File(f2).NT(Cand(1));
        NT4 = File(f2).NT(Cand(2));

        [d,ang,t1,t2,cp] = zIsoDiscrepancy(NT1,NT2,NT3,NT4);
        ID(k,m)  = d;
        Ang(k,m) = ang^2;
        T1(k,m)  = t1*t1';
        T2(k,m)  = t2*t2';
        CP(k,m)  = cp^2;
      end
    end

    PD = sqrt(PD)/2;            % finish discrepancy calculation

  if 0 > 1,
    bigm = max(max(PD));
    if bigm > 1,
      fprintf('%6.2f maximum pair discrepancy\n',bigm);
    else
      fprintf('\n');
    end
  end

  if Verbose > 0,
    fprintf(' %4d coplanar. ', sum(diag(PD) == 0));
  end

  PD = PD + PD';
  ID = ID + ID';

  rs = sum(PD);                              % total geometric discrepancy to all other instances
  [y,i] = sort(rs);                          % sort by centrality

  qs = sum(ID);                              % total isodiscrepancy to all other instances
  [z,k] = sort(qs);                          % sort by centrality

  f = List(i(1),3);

  if Verbose > 0,
    fprintf('Coplanar value of exemplar is %7.4f\n', full(File(f).Coplanar(List(i(1),1),List(i(1),2))));
  end
