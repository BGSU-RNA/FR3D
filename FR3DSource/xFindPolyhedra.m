% xFindHexahedra.m returns a list of indices (i,j,k,m) for which 
% A(i,j) > 0, B(k,i) > 0, C(k,j) > 0, D(m,i) > 0, E(m,j) > 0, F(m,k) > 0
% G(r,i) > 0, H(r,j) > 0, I(r,m) > 0, J(r,k) > 0

function [TList,SS] = xFindPolyhedra(Model,num,PS)

SS = [];
A=PS{2,1};
if num>2
    Cutoff=Model.SSCutoff;
    B=PS{3,1};
    C=PS{3,2};
    if num>3
        D=PS{4,1};
        E=PS{4,2};
        F=PS{4,3};
        if num>4
            G=PS{5,1};
            H=PS{5,2};
            I=PS{5,3};
            J=PS{5,4};
            if num>5
                K=PS{6,1};
                L=PS{6,2};
                M=PS{6,3};
                N=PS{6,4};
                O=PS{6,5};
                if num>6
                    P=PS{7,1};
                    Q=PS{7,2};
                    R=PS{7,3};
                    S=PS{7,4};
                    T=PS{7,5};
                    U=PS{7,6};
                    if num>7
                        V=PS{8,1};
                        W=PS{8,2};
                        X=PS{8,3};
                        Y=PS{8,4};
                        Z=PS{8,5};
                        AA=PS{8,6};
                        BB=PS{8,7};
                        if num>8
                            CC=PS{9,1};
                            DD=PS{9,2};
                            EE=PS{9,3};
                            FF=PS{9,4};
                            GG=PS{9,5};
                            HH=PS{9,6};
                            II=PS{9,7};
                            JJ=PS{9,8};
                        end
                    end
                end
            end
        end
    end
end

[j,i] = find(A);
LL = 64000;                                % initial size for list storage
TList = uint16(zeros(LL,num));
count = 0;
switch num
    case 2
        TList   = [i j];
        count=length(i);
    case 3
        [TList,count]=Case3(Cutoff,TList,LL,i,j,count,A,B,C); 
        
    case 4
        [TList,count]=Case4(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F); 

    case 5
        [TList,count]=Case5(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J); 
        
    case 6
        [TList,count]=Case6(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O);

    case 7
        [TList,count]=Case7(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U);
        
    case 8
        [TList,count]=Case8(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB); 
        
    case 9
        SS = zeros(LL,1);
        [TList,count,SS]=Case9(Cutoff,TList,SS,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ);
        SS = SS(1:count,1);
end
TList = TList(1:count,:);

%---------------------------------------------------------------
function [TList,count]=Case3(Cutoff,TList,LL,i,j,count,A,B,C)
    

TList = uint16(zeros(LL,3));
count = 0;
for n=1:length(i),
  in = i(n);
  jn = j(n);  
  k = find(B(:,in) .* C(:,jn));
  K = length(k);
  if count + K > LL - 10000,                      % memory management
    AddLL = min(LL,1000000);
    TList = [TList; uint16(zeros(AddLL,3))];
    LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
  end
  for p = 1:K,
    kp=k(p);
    SSks= A(jn,in) + B(kp,in) + C(kp,jn);
    count = count + 1;
    TList(count,:) = [in jn kp];
  end
end


%---------------------------------------------------------------
function [TList,count]=Case4(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F)

  for n=1:length(i),
      in = i(n);
      jn = j(n);
      k = find(B(:,in) .* C(:,jn));
      m = find(D(:,in) .* E(:,jn));
      [mm,kk]  = find(F(m,k));
      if count + length(kk) > LL - 10000,                      % memory management
        AddLL = min(LL,1000000);
        TList = [TList; uint16(zeros(AddLL,4))];
        LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
      end
      for p = 1:length(kk),
          kkkp=k(kk(p));
          mmmp=m(mm(p));
        SSks= A(jn,in) + B(kkkp,in) + C(kkkp,jn) ...
            + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
        if SSks<Cutoff(4)
            count = count + 1;
            TList(count,:) = [in jn kkkp mmmp];
        end
      end
  end



  %---------------------------------------------------------------
function [TList,count]=Case5(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J)

  for n=1:length(i),
      in = i(n);
      jn = j(n);      
      k = find(B(:,in) .* C(:,jn));
      m = find(D(:,in) .* E(:,jn));
      [mm,kk]  = find(F(m,k));
      r = find(G(:,in) .* H(:,jn));
      for p = 1:length(kk),
            kkkp = k(kk(p));
            mmmp = m(mm(p));          
           SSkm= A(jn,in) + B(kkkp,in) + C(kkkp,jn) ...
            + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
            if SSkm<Cutoff(4)
            rr = find(I(r,kkkp) .* J(r,mmmp));
                for q = 1:length(rr),
                    rrrq=r(rr(q));
                    SSrs = SSkm + G(rrrq,in)+H(rrrq,jn)+I(rrrq,kkkp)+J(rrrq,mmmp);
                    if SSrs<Cutoff(5)
                        if count + length(kk) > LL - 10000,                      % memory management
                          AddLL = min(LL,1000000);
                          TList = [TList; uint16(zeros(AddLL,5))];
                          LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
                        end
                        count = count + 1;
                        TList(count,:) = [in jn kkkp mmmp rrrq];
                    end
                end
            end
      end
  end % end for n
%---------------------------------------------------------------
function [TList,count]=Case6(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)

    for n=1:length(i),
        in = i(n);
        jn = j(n);
        SSij = A(jn,in);
        k = find(B(:,in) .* C(:,jn));
        m = find(D(:,in) .* E(:,jn));
        [mm,kk] = find(F(m,k));
        if length(mm) > 0,
            r = find(G(:,in) .* H(:,jn));
            s = find(K(:,in) .* L(:,jn));

            for p = 1:length(kk),
            kkkp = k(kk(p));
            mmmp = m(mm(p));
            SSkm = SSij + B(kkkp,in) + C(kkkp,jn) ...
            + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
                if SSkm < Cutoff(4)
                    rr = find(I(r,kkkp) .* J(r,mmmp));
                    ss = find(M(s,kkkp) .* N(s,mmmp));
                    [sr,rs] = find(O(s(ss),r(rr)));
                    for q = 1:length(rs),
                        rrrrsq=r(rr(rs(q)));
                        ssssrq=s(ss(sr(q)));
                        if count + length(rs) > LL - 10000,                      % memory management
                          AddLL = min(LL,1000000);
                          TList = [TList; uint16(zeros(AddLL,6))];
                          LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
                        end                   
                        SSrs = SSkm + G(rrrrsq,in)+H(rrrrsq,jn)+I(rrrrsq,kkkp)+J(rrrrsq,mmmp)+...
                            K(ssssrq,in)+L(ssssrq,jn)+M(ssssrq,kkkp)+N(ssssrq,mmmp)+O(ssssrq,rrrrsq);
                        if SSrs < Cutoff(6)
                            count = count + 1;
                            TList(count,:) = [in jn kkkp mmmp rrrrsq ssssrq];
                        end
                    end
                end
            end
        end
    end
  %---------------------------------------------------------------
function [TList,count]=Case7(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U)
    for n=1:length(i),
      in = i(n);
      jn = j(n);
      k = find(B(:,in) .* C(:,jn));
      m = find(D(:,in) .* E(:,jn));
      [mm,kk]  = find(F(m,k));
        r = find(G(:,in) .* H(:,jn));
        s = find(K(:,in) .* L(:,jn));
        t = find(P(:,in) .* Q(:,jn));
        for p = 1:length(kk),    
            kkkp = k(kk(p));
            mmmp = m(mm(p));            
            SSkm= A(jn,in) + B(kkkp,in) + C(kkkp,jn) ...
                + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
            if SSkm<Cutoff(4)          
                rr = find(I(r,kkkp) .* J(r,mmmp));
                ss = find(M(s,kkkp) .* N(s,mmmp));
                tt = find(R(t,kkkp) .* S(t,mmmp));
                [sr,rs] = find(O(s(ss),r(rr)));
                 for q=1:length(rs)   
                    rrrrsq=r(rr(rs(q)));
                    ssssrq=s(ss(sr(q)));                     
                    SSrs = SSkm + G(rrrrsq,in)+H(rrrrsq,jn)+I(rrrrsq,kkkp)+J(rrrrsq,mmmp)+...
                        K(ssssrq,in)+L(ssssrq,jn)+M(ssssrq,kkkp)+N(ssssrq,mmmp)+O(ssssrq,rrrrsq);
                    if SSrs < Cutoff(6) 
                        trs = find(T(t(tt),rrrrsq) .* U(t(tt),ssssrq));    
                        for e=1:length(trs)
                            ttttrse=t(tt(trs(e)));
                            SSts=SSrs+P(ttttrse,in)+Q(ttttrse,jn)+R(ttttrse,kkkp)+S(ttttrse,mmmp)+T(ttttrse,rrrrsq)+U(ttttrse,ssssrq);
                            if SSts<Cutoff(7)
                                if count + length(trs) > LL - 10000,                      % memory management
                                  AddLL = min(LL,1000000);
                                  TList = [TList; uint16(zeros(AddLL,7))];
                                  LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
                                end        
                                count = count + 1;
                                TList(count,:) = [in jn kkkp mmmp rrrrsq ssssrq t(tt(trs(e)))];
                            end
                        end
                    end
                 end
            end
        end
    end
  %---------------------------------------------------------------
  function [TList,count]=Case8(Cutoff,TList,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB)
    for n=1:length(i),
      in = i(n);
      jn = j(n);        
      k = find(B(:,in) .* C(:,jn));
      m = find(D(:,in) .* E(:,jn));
      [mm,kk]  = find(F(m,k));
        r = find(G(:,in) .* H(:,jn));
        s = find(K(:,in) .* L(:,jn));
        t = find(P(:,in) .* Q(:,jn));
        h = find(V(:,in) .* W(:,jn));
        if count + length(kk) > length(TList(:,1)),                  % memory management
            TList = [TList; uint16(zeros(LL,8))];
            LL = min(2*LL,1000000);
        end
        for p = 1:length(kk),
            kkkp = k(kk(p));
            mmmp = m(mm(p));            
            SSkm= A(jn,in) + B(kkkp,in) + C(kkkp,jn) ...
                + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
            if SSkm<Cutoff(4)            
            
                rr = find(I(r,kkkp) .* J(r,mmmp));
                ss = find(M(s,kkkp) .* N(s,mmmp));
                tt = find(R(t,kkkp) .* S(t,mmmp));
                hh = find(X(h,kkkp) .* Y(h,mmmp));
                [sr,rs] = find(O(s(ss),r(rr)));

                 for q=1:length(rs)
                    rrrrsq=r(rr(rs(q)));
                    ssssrq=s(ss(sr(q)));                     
                    SSrs = SSkm + G(rrrrsq,in)+H(rrrrsq,jn)+I(rrrrsq,kkkp)+J(rrrrsq,mmmp)+...
                        K(ssssrq,in)+L(ssssrq,jn)+M(ssssrq,kkkp)+N(ssssrq,mmmp)+O(ssssrq,rrrrsq);
                    if SSrs < Cutoff(6) 

                        trs = find(T(t(tt),rrrrsq) .* U(t(tt),ssssrq));    
                        hrs = find(Z(h(hh),rrrrsq) .* AA(h(hh),ssssrq));
                        [ht,th] = find(BB(h(hh(hrs)),t(tt(trs))));
                        for e=1:length(ht)
                            ttttrsthe=t(tt(trs(th(e))));
                            hhhhrshte=h(hh(hrs(ht(e))));
                            SSts=SSrs+P(ttttrsthe,in)+Q(ttttrsthe,jn)+R(ttttrsthe,kkkp)+S(ttttrsthe,mmmp)+T(ttttrsthe,rrrrsq)+U(ttttrsthe,ssssrq)+...
                                V(hhhhrshte,in)+W(hhhhrshte,jn)+X(hhhhrshte,kkkp)+Y(hhhhrshte,mmmp)+Z(hhhhrshte,rrrrsq)+AA(hhhhrshte,ssssrq)+BB(hhhhrshte,ttttrsthe);
                            if SSts<Cutoff(8)
                                if count + length(ht) > LL - 10000,                      % memory management
                                  AddLL = min(LL,1000000);
                                  TList = [TList; uint16(zeros(AddLL,8))];
                                  LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
                                end
                                count = count + 1;
                                TList(count,:) = [in jn kkkp mmmp rrrrsq ssssrq ttttrsthe hhhhrshte];
                            end
                        end
                    end
                 end
            end
        end
    end
  %---------------------------------------------------------------
    function [TList,count,SS]=Case9(Cutoff,TList,SS,LL,i,j,count,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ)

    for n=1:length(i),                         % loop through pairs
      in = i(n);                               % indices of this pair
      jn = j(n);        
      k = find(B(:,in) .* C(:,jn));            % k that agree with pair
      m = find(D(:,in) .* E(:,jn));            % m that agree with pair
      [mm,kk]  = find(F(m,k));                 % k and m that match too
        r = find(G(:,in) .* H(:,jn));
        s = find(K(:,in) .* L(:,jn));
        t = find(P(:,in) .* Q(:,jn));
        h = find(V(:,in) .* W(:,jn));
        u = find(CC(:,in) .* DD(:,jn));
        for p = 1:length(kk),                  % loop through tetrahedra
            kkkp = k(kk(p));                   % tetrahedron in jn kkkp mmmp
            mmmp = m(mm(p));            
            SSkm = A(jn,in) + B(kkkp,in) + C(kkkp,jn) ...
                 + D(mmmp,in) + E(mmmp,jn) + F(mmmp,kkkp);
            if SSkm<Cutoff(4)            
                rr = find(I(r,kkkp) .* J(r,mmmp));
                ss = find(M(s,kkkp) .* N(s,mmmp));
                tt = find(R(t,kkkp) .* S(t,mmmp));
                hh = find(X(h,kkkp) .* Y(h,mmmp));
                uu = find(EE(u,kkkp).* FF(u,mmmp));
                [sr,rs] = find(O(s(ss),r(rr)));
                 for q=1:length(rs)
                    rrrrsq=r(rr(rs(q)));
                    ssssrq=s(ss(sr(q)));                     
                    SSrs = SSkm + G(rrrrsq,in)+H(rrrrsq,jn)+I(rrrrsq,kkkp)+J(rrrrsq,mmmp)+...
                        K(ssssrq,in)+L(ssssrq,jn)+M(ssssrq,kkkp)+N(ssssrq,mmmp)+O(ssssrq,rrrrsq);
                    if SSrs < Cutoff(6) 
                        trs = find(T(t(tt),rrrrsq) .* U(t(tt),ssssrq));    
                        hrs = find(Z(h(hh),rrrrsq) .* AA(h(hh),ssssrq));
                        urs = find(GG(u(uu),rrrrsq) .* HH(u(uu),ssssrq));
                        [ht,th] = find(BB(h(hh(hrs)),t(tt(trs))));
                        for e=1:length(ht)
                            ttttrsthe=t(tt(trs(th(e)))); 
                            hhhhrshte=h(hh(hrs(ht(e))));
                               ust = find(II(u(uu(urs)),ttttrsthe) .* JJ(u(uu(urs)),hhhhrshte));                          
                               for ee=1:length(ust)
                                   uuuursustee=u(uu(urs(ust(ee))));
                                    SSts=SSrs+P(ttttrsthe,in)+Q(ttttrsthe,jn)+R(ttttrsthe,kkkp)+S(ttttrsthe,mmmp)+T(ttttrsthe,rrrrsq)+U(ttttrsthe,ssssrq)+...
                                        V(hhhhrshte,in)+W(hhhhrshte,jn)+X(hhhhrshte,kkkp)+Y(hhhhrshte,mmmp)+Z(hhhhrshte,rrrrsq)+AA(hhhhrshte,ssssrq)+BB(hhhhrshte,ttttrsthe)+...
                                        CC(uuuursustee,in)+DD(uuuursustee,jn)+EE(uuuursustee,kkkp)+FF(uuuursustee,mmmp)+GG(uuuursustee,rrrrsq)+...
                                        HH(uuuursustee,ssssrq)+II(uuuursustee,ttttrsthe)+JJ(uuuursustee,hhhhrshte);
                                    if SSts < Cutoff(9)
                                        if count + length(ust) > LL - 10000,                      % memory management
                                          AddLL = min(LL,1000000);
                                          TList = [TList; uint16(zeros(AddLL,9))]; 
                                          SS    = [SS; zeros(AddLL,1)];            
                                          LL    = length(TList(:,1));
if count > 1000000,
  fprintf('Found %2d million candidates so far\n',fix(count/1000000));
end
                                        end
                                        count = count + 1;
                                        TList(count,:) = [in jn kkkp mmmp rrrrsq ssssrq ttttrsthe hhhhrshte uuuursustee];
                                        SS(count,1)    = SSts;
                                    end
                               end
                        end
                    end
                 end
            end
        end
    end
  %---------------------------------------------------------------

