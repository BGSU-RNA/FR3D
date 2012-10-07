%  xFindPermutation(b);

function [Perm] = xFindPermutation(NNZ);

[m,n] = size(NNZ);

count = 0;
for i=2:min(9,n)
  for j=1:(i-1)
    count=count+1;
    LEN(count,1)=NNZ(i,j);
    LEN(count,2:3)=[i,j];
  end
end

Loc=1;
i=1:n;
while Loc+1<n 
    [a,b]=sortrows(LEN,1);
    [LEN,i]=MakeMinFirst(LEN,b,i,Loc);
    Loc=Loc+2;
end

Perm=i;

%---------------------------------------------------
function [LEN,i]=MakeMinFirst(LEN,b,i,Loc)

% switch first number
k=find(i==LEN(b(1),2));
t=i(Loc);
i(Loc)=LEN(b(1),2); 
i(k)=t;

% switch 2nd number
k=find(i==LEN(b(1),3));
t=i(Loc+1);
i(Loc+1)=LEN(b(1),3); 
i(k)=t;

c=find(LEN(:,2)~=LEN(b(1),3) & LEN(:,2)~=LEN(b(1),2) & LEN(:,3)~=LEN(b(1),3) & LEN(:,3)~=LEN(b(1),2));
LEN=LEN(c,:);





