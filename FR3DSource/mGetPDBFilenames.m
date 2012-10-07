%mGetPDBfilenames %By Ali Mokdad - April 14 2006
%This looks in the folder <PDBFiles> for files with extensions .pdb-1repetitions

% Extended to search Matlab's path for PDB files and to use the native
% file separator rather than \  - CLZ 2006-07-18

function [s,snolist] = mGetPDBFilenames(void)

a = dir(['PDBFiles' filesep '*.pdb']);
a = [a; dir(['PDBFiles' filesep '*.PDB'])];         % for the Mac
a = [a; dir(['PrecomputedData' filesep '*.mat'])];
a = [a; dir(['PrecomputedData' filesep '*.MAT'])];  % for the Mac

p = path;                                           % search Matlab's path
c = [0 strfind(p,pathsep) length(p)+1];
for i=1:length(c)-1,
  a = [a; dir([p(c(i)+1:c(i+1)-1) filesep '*.pdb'])];
end

if ~isempty(a)
    for i=1:length(a)
      temp{i} = regexprep(a(i).name,'.mat|.MAT|.PDB1|.pdb1|.PDB2|.pdb2|.PDB3|.pdb3|.PDB4|.pdb4|.PDB5|.pdb5|.PDB6|.pdb6|.PDB7|.pdb7|.PDB8|.pdb8','');   % strip extensions and spaces
      temp{i} = regexprep(temp{i},'.pdb|.PDB ','');   % strip extensions and spaces
      temp{i} = regexprep(temp{i},'_small','');     % don't list small's too
      if ~isempty(strfind(temp{i},'_list')),
        temp{i} = [' ' temp{i}];                    % lists appear first
      end
    end

    temp = [' AllFiles_list' temp];
    
    [t,i] = sort(upper(temp));                       % sort list of PDB names
    temp = temp(i);

    for i=1:length(temp),
      temp{i} =regexprep(temp{i},' ','');           % strip spaces from lists
    end

    s{1,1}=temp{1};
    count=2;
    for i=2:length(temp)
      if ~strcmpi(temp{i},s{count-1,1})              % remove duplicates
        s{count,1}=temp{i};
        count=count+1;
      end
    end

    count = 1;

    for i = 1:length(s(:,1)),
      if isempty(strfind(s{i,1},'_list')),
        snolist{count,1} = s{i,1};
        count = count + 1;
      end
    end

else
  s = [];
  snolist = []; 
end

