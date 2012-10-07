% xAddNucleotide adds a nucleotide according to the screen matrices S

function [NewList,NewSS] = xAddNucleotide(Model,List,S,SS,r)

L = 100000;

[s,t] = size(List);

NewList = uint16(zeros(L,t+1));
NewSS   = zeros(L,1);

A = S{r,1};                % for some reason, this is MUCH faster
B = S{r,2};
C = S{r,3};
D = S{r,4};

if r >  5,  E = S{r,5}; end
if r >  6,  F = S{r,6}; end
if r >  7,  G = S{r,7}; end
if r >  8,  H = S{r,8}; end
if r >  9,  I = S{r,9}; end

c = 0;                                     % number of candidates

for n=1:s,
  a = List(n,1:(r-1));                     % a is the ith candidate
  switch r                                 % search for one more nucleotide
    case 5,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4));
    case 6,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5));
    case 7,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5)) .*  ...
          F(:,a(6));
    case 8,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5)) .*  ...
          F(:,a(6)) .* G(:,a(7));
    case 9,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5)) .*  ...
          F(:,a(6)) .* G(:,a(7)) .* H(:,a(8));
    case 10,
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5)) .*  ...
          F(:,a(6)) .* G(:,a(7)) .* H(:,a(8)) .* I(:,a(9));
    otherwise
      h = A(:,a(1)) .* B(:,a(2)) .* C(:,a(3)) .* D(:,a(4)) .* E(:,a(5)) .*  ...
          F(:,a(6)) .* G(:,a(7)) .* H(:,a(8)) .* I(:,a(9));
      for q = 10:(r-1), 
        h = h .* S{r,q}(:,a(q));           % slower than hard-coding
      end
  end

  m = find(h);                             % indices of new nucleotides

  if c + length(m) > L,                    % memory management; faster
    NewList = [NewList; uint16(zeros(L,t+1))];
    NewSS   = [NewSS; zeros(L,1)];
    L = 2*L;
  end

  for i=1:length(m),                       % append these new candidates
    c = c + 1;
    NewList(c,:) = [a m(i)];
    NewSS(c,1)   = SS(n,1);
  end
end

NewList = NewList(1:c,:);
NewSS   = NewSS(1:c,:);

if Model.Geometric > 0,
  Acceptable = uint16(zeros(size(NewList(:,1))));
  for n = 1:length(NewList(:,1)),
    a = NewList(n,:);
    b = a(r);
    y = NewSS(n,1) + A(b,a(1)) + B(b,a(2)) + C(b,a(3)) + D(b,a(4));
    switch r,
      case  6, y = y + E(b,a(5));
      case  7, y = y + E(b,a(5)) + F(b,a(6));
      case  8, y = y + E(b,a(5)) + F(b,a(6)) + G(b,a(7));
      case  9, y = y + E(b,a(5)) + F(b,a(6)) + G(b,a(7)) + H(b,a(8));
      case 10, y = y + E(b,a(5)) + F(b,a(6)) + G(b,a(7)) + H(b,a(8)) + ...
                       I(b,a(9));
      case 11, 
        y = y + E(b,a(5)) + F(b,a(6)) + G(b,a(7)) + H(b,a(8)) + I(b,a(9));
        for q = 10:(r-1),
          y = y + S{r,q}(b,a(q));
        end
    end

    if y < Model.SSCutoff(r),
      Acceptable(n) = 1;
    end

  end

  OK      = find(Acceptable);  
  NewList = NewList(OK,:);
  NewSS   = NewSS(OK,:);

%  fprintf('Retained %6d, rejected %5d candidates on nucleotide %1d\n', length(OK), length(Acceptable)-length(OK), r);
end

