% zSelectPairs(File,Param) returns a list of indices of selected bases and

% Param(m,1) = interaction code for mth pair
% Param(m,2) = paircode

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

function [List] = zSelectPairs(File,Param,Decimal)

[s,t] = size(Param);

if t > 2,
  for m = 1:length(Param(:,1)),
    pc = 4*(Param(m,3)-1) + Param(m,2);
    Param(m,4) = Param(m,3);
    Param(m,3) = Param(m,2);
    Param(m,2) = pc;
  end
else
  for m = 1:length(Param(:,1)),
    pc = Param(m,2);
    Param(m,3) = mod(pc-1,4)+1;
    Param(m,4) = fix((pc-1)/4)+1;
  end
end

for m = 1:length(Param(:,1)),
  if (Param(m,2) == 7) || ((Param(m,2) == 10) && (abs(Param(m,1)) < 14)),
    Param(m,1) = -Param(m,1);
  end
end

List = [];

for f = 1:length(File),                   % Loop through each file
 Codes = cat(1,File(f).NT(:).Code);

 if Decimal == 1,                         % Use decimal places in category 
   G = File(f).Edge;
 else                                     % Ignore decimal places in category
   G = fix(File(f).Edge);
 end

 N = length(File(f).NT);                     % Number of nucleotides in File

 for m = 1:length(Param(:,1)),
   [i,j] = find(G == Param(m,1));

   k = find((Codes(i) == Param(m,3)) .* (Codes(j) == Param(m,4)));

   List = [List; [i(k) j(k) f*ones(length(k),1)]];
 end

end
