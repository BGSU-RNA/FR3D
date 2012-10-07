% FR3D_InteractionMap(File) puts a dot at the center of each nucleotide
% in File, colored by the number of interactions the nucleotide has (blue
% means zero, red is the most).  The C1' atom of each nucleotide is
% connected to the C1' atom of the next nucleotide by a black line.
% Interactions are shown between centers of nucleotides.  WC-WC
% interactions (category 1) are shown by red lines, non-canonical planar
% interactions by green, typical sequential stacking (category 21) by
% dark blue, and typical cross-strand stacking (categories 22, 23) by light
% blue.
% File can be, for example, '1s72' or File(3).
% The vector DL, if specified, controls what is displayed:
%  DL(1) = 1;                                 % plot geometric center
%  DL(2) = 1;                                 % show cWW bonds
%  DL(3) = 1;                                 % show other planar bonds
%  DL(4) = 1;                                 % show 35 and 53 stacking
%  DL(5) = 1;                                 % show 33 and 55 stacking
%  DL(6) = 1;                                 % connect glycosidic atoms

function [void] = FR3D_InteractionMap(File,DL,NTList)

if nargin < 1
  fprintf('Please specify a filename, e.g. FR3D_InteractionMap(''1s72'')\n');
end

% if File is a text string (filename), load the file and display

if strcmp(class(File),'char'),
  Filename = File;
  File = zGetNTData(Filename,0);
end

% if NTList is a cell array of numbers, look up the indices

if nargin > 2,

  if strcmp(class(NTList),'char'),
    NTList = {NTList};
  end

  if strcmp(class(NTList),'cell'),
    Indices = zIndexLookup(File,NTList);
  else
    Indices = NTList;
  end

else

  Indices = 1:length(File.NT);

end

if nargin < 2
  DL(1) = 1;                                 % plot geometric center
  DL(2) = 1;                                 % show cWW bonds
  DL(3) = 1;                                 % show other planar bonds
  DL(4) = 1;                                 % show 35 and 53 stacking
  DL(5) = 1;                                 % show 33 and 55 stacking
  DL(6) = 1;                                 % connect glycosidic atoms
end

%-----------------------------------------------------------------------
figure(1)
clf

set(gcf,'Renderer','OpenGL');

for a = 1:length(Indices),                          % Loop through nucleotides
  k = Indices(a);
  i  = abs(File.Edge(k,Indices));
  ii = find((i < 30) .* (i > 0));
  iii = Indices(ii);
  ni(k) = length(iii);                           % number of interactions

  e = File.NT(k).Center;
  if DL(1) > 0,
    scatter3(e(1),e(2),e(3),6,ni(k)+1,'filled');     % plot geometric center
  else
    scatter3(e(1),e(2),e(3),1,ni(k)+1,'filled');     % plot dots at least
  end
  hold on

  c = 0;

  for j = 1:length(iii),
    if iii(j) < k,                            % only plot once
      f = File.NT(iii(j)).Center;
      inter = abs(fix(File.Edge(k,iii(j))));
      if inter == 1,
        if DL(2) > 0,
          plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'b');  % WC-WC
        end
        c = c + 1;
      elseif abs(inter) < 15,
        if DL(3) > 0,
          plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'c');  % other planar
        end
        c = c + 1;
      elseif (inter == 21),
        if DL(4) > 0,
          plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'y');  % 35 or 53 stacking
        end
        c = c + 1;
      elseif (inter == 22) | (inter == 23),
        if DL(5) > 0,
          plot3([e(1) f(1)],[e(2) f(2)],[e(3) f(3)],'y');  % 33 or 55 stacking
        end
        c = c + 1;
      end
    end
  end

  if DL(6) > 0,
    g = File.NT(k).Fit(1,:);                         
    h = File.NT(max(1,k-1)).Fit(1,:);
    if sum((g-h).^2) < 200,                           % reasonable distance
      plot3([g(1) h(1)],[g(2) h(2)],[g(3) h(3)],'k')  % connect glycosidics
    end
  end

end

set(gcf,'Renderer','painters');
axis off

rotate3d on