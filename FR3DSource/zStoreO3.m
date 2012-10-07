
function [File] = zStoreO3(File)

for f = 1:length(File),
  File(f).NT(1).Sugar(13,:) = [Inf Inf Inf];
  for i = 2:length(File(f).NT),
    if File(f).NT(i).Chain == File(f).NT(i-1).Chain,
      File(f).NT(i).Sugar(13,:) = File(f).NT(i-1).Sugar(5,:);
    else
      File(f).NT(i).Sugar(13,:) = File(f).NT(i).Sugar(10,:); % use phosphorus
    end
  end
end
