% FR3D_interaction_list produces lists of all FR3D annotated interactions

% To run from the command line, matlab -r FR3D_interaction_list('1s72')
% To run from the command line, matlab -r FR3D_interaction_list('1s72','C:\Users\zirbel\Documents\FR3D')

% 1S72_interactions_FR3D.txt contains all interactions except near
% 1S72_near_interactions_FR3D.txt contains near interactions
% 1S72_stacking_FR3D.txt contains stacking interactions
% 1S72_basepairs_FR3D.txt contains basepairing interactions
% 1S72_base_backbone_FR3D.txt contains base-phosphate and base-ribose

% Note:  basepairs also include bifurcated (bif), water-inserted (wat),
%        ribose basepairs (rib) and perpendicular interactions (perp)
% Only the bifurcated pairs have been written about.

% Note:  base-phosphate interactions were discussed in a published article.
%        base-ribose interactions follow the same classification routine,
%        but use the O2', O3', and O4' atoms on the ribose sugar

function [void] = zAnalyzedFilesHTML(File,datapath)

path(path,[pwd filesep 'FR3DSource']);

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char') || strcmp(class(File),'cell'),
  Filename = File;
  File = zAddNTData(Filename,0,[],1);
end

for f = 1:length(File),
 if length(File(f).NT) > 1,
  clear chainoffset

  FN = upper(File(f).Filename);
  PDBFN = File(f).PDBFilename;

  Vers = num2str(File(f).ClassVersion);

  DataHeader1 = sprintf('# PDB_ID_FR3D_Version_%s',Vers);

  DataHeader2 = sprintf('PDB_ID\tInteraction\tNucleotide_1_Base\tNucleotide_1_PDB_Number\tNucleotide_1_Chain\tNucleotide_1_Sequence_Position\tNucleotide_2_Base\tNucleotide_2_PDB_Number\tNucleotide_2_Chain\tNucleotide_2_Sequence_Position');

  % -------------------------------------------------- Set paths

  warning off

  if nargin == 1,
    datapath = [pwd filesep];
  end

  if ~exist([datapath FN],'dir'), 
    mkdir([datapath FN]); 
  end

  % ------------------------------------------- Write chains and sequences

  Chain = cat(2,File(f).NT.Chain);
  U     = unique(Chain);
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    chainoffset(u) = min(i);                    % first index in this chain
  end

  % ------------------------------------------- Produce interaction list

  c = 1;                                    % counter for interactions

  IText{1} = '';
  DText{1} = '';
  InterType = [];

  E   = File(f).Edge;
  BPh = File(f).BasePhosphate;
  BR  = File(f).BaseRibose;
  BC  = File(f).Covalent;                    % covalent connections

  for i = 1:File(f).NumNT,
    N1 = File(f).NT(i);

    % ------------------------------------- Find basepairing, stacking
    j = find(E(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));
      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code),r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;         % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;      % position of j(k) in chain u

      T = zEdgeText(File(f).Edge(i,j(k)),0,N1.Code,N2.Code);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      InterType(c) = abs(File(f).Edge(i,j(k)));

      c = c + 1;
    end

    % ------------------------------------- Find base phosphate interactions
    j = find(BPh(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBasePhosphateText(BPh(i,j(k)),1), r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zBasePhosphateText(BPh(i,j(k)),1);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      if abs(BPh(i,j(k))) > 100,
        InterType(c) = abs(BPh(i,j(k)));      % near interaction
      elseif i == j(k),
        InterType(c) = 200.1;                 % self interaction
      else
        InterType(c) = 200;                   % non-self interaction
      end

      c = c + 1;
    end

    % ------------------------------------- add base ribose interactions
    j = find(BR(i,:));
    for k = 1:length(j),
      
      N2 = File(f).NT(j(k));
      r = sprintf('%4d', full(File(f).Crossing(i,j(k))));

      IText{c} = sprintf('%s%4s(%s) - %s%4s(%s) - %7s - %s', N1.Base, N1.Number, N1.Chain, N2.Base, N2.Number, N2.Chain, zBaseRiboseText(BR(i,j(k)),1), r);

      u  = find(U==N1.Chain);              % which chain i is in
      ii = i - chainoffset(u) + 1;        % position of i in chain u
      u  = find(U==N2.Chain);              % which chain j(k) is in
      jj = j(k) - chainoffset(u) + 1;     % position of j(k) in chain u

      T = zBaseRiboseText(BR(i,j(k)),1);

      DText{c} = sprintf('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n', File(f).Filename, T, N1.Base, N1.Number, N1.Chain, ii, N2.Base, N2.Number, N2.Chain, jj);

      if abs(BR(i,j(k))) > 100,
        InterType(c) = abs(BR(i,j(k)));      % near interaction
      elseif i == j(k),
        InterType(c) = 200.1;                 % self interaction
      else
        InterType(c) = 200;                   % non-self interaction
      end

      c = c + 1;
    end
  end

  % --------------------------------------------- Write FN_interactions file  

  fid = fopen([datapath FN '_interactions_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find((InterType < 30) + (InterType >= 200));  % exclude near pairs
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

  % ----------------------------------------- Write FN_near_interactions file  


  fid = fopen([datapath FN '_near_interactions_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find((InterType > 100) .* (InterType < 200) .* (fix(InterType) ~= 114));  % exclude near pairs, BPh
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

  % ----------------------------------------------- Write FN_basepairs file  


  fid = fopen([datapath FN '_basepairs_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find(InterType < 14);                    % exclude Rib pairs
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% ----------------------------------------------- Write FN_stacking file  

  fid = fopen([datapath FN '_stacking_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find((InterType > 19) .* (InterType < 25));
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

% --------------------------------------------- Write FN_base_backbone file  

  fid = fopen([datapath FN '_base_backbone_FR3D.txt'],'w'); % open for writing
  fprintf(fid,'%s\n',DataHeader1);
  fprintf(fid,'# Chain\tNucleotide_sequence_in_chain\n');
  for u = 1:length(U),
    i = find(Chain == U(u));                    % NTs in chain U(u)
    fprintf(fid,'# %s\t%s\n',U(u),cat(2,File(f).NT(i).Base));
  end
  fprintf(fid,'%s\n',DataHeader2);
  k = find(fix(InterType) == 200);               % all BPh and BR interactions
  for i = 1:length(k),
    fprintf(fid,'%s',strrep(DText{k(i)},' ',''));
  end
  fclose(fid);

  clear IText HText DText

 end
end


