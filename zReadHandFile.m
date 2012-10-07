% zReadHandFile reads in previously created expert classifications.
% These are stored in the following fields of File:

% HandClass    m x 1 array of expert classification categories
% Comment      m row cell array of comments on each pair
% CI           NumNT x NumNT matrix of indices to the previous two arrays

% Thus, if you want to know what the expert classification is of the
% interaction between bases with indices p and q, you look at
% Comment{CI(p,q)} and HandClass(CI(p,q))

function [File] = zReadHandFile(File)

HandClass = [];                         % hand category read from file
Comment  = [];                         % comments read from file
CI        = sparse(File.NumNT,File.NumNT);        % matrix of indices

if exist(strcat(File.Filename,'.hand'),'file') > 0,
  [r_A, r_B, r_H, r_C, r_D, r_I, r_E, r_F, r_G] ...
  = textread(strcat(File.Filename,'.hand'),'%2s %5s %2s %2s %5s %2s %6.2f %6.2f %[^\n]','endofline','\n');

  NuclNum1   = r_B;
  NuclNum2   = r_D;
  HandClass  = r_F;
  Comment    = r_G;

  Numbers   = cat(1,{File.NT(:).Number});

  for a=1:length(HandClass),
    p = find(ismember(Numbers,NuclNum1{a}));  % get index of first nucleotide
    q = find(ismember(Numbers,NuclNum2{a}));  % get index of second nucleotide

% in the future, check for p and q having length > 1
% also check for chain number here

    if (p>0) & (q>0),
      if CI(p,q) > 0,
        fprintf('Nucleotides %5s and %5s ',NuclNum1{a}, NuclNum2{a});
        fprintf('are hand classified twice.\n');
      end
      CI(p,q) = a;                        % points pair to row of Comment
      CI(q,p) = a;                        % points pair to row of Comment
      if (File.Distance(p,q) == 0),
        fprintf('Nucleotides %5s and %5s ',NuclNum1{a}, NuclNum2{a});
        fprintf('are more than 15 Angstroms apart.\n');
      end
    else
      fprintf('Could not find one of nucleotides ');
      fprintf('%s and %s in class file.\n', NuclNum1{a}, NuclNum2{a});
    end
  end
  fprintf('Read  %s\n', [File.Filename '.hand']);
end

% Add fields to File ---------------------------------------------------------

File.HandClass = HandClass;
File.Comment   = Comment;
File.CI        = CI;
