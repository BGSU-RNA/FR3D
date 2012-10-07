% zBaseRiboseInteractions adapts zPhosphateInteractions to detect
% base-ribose interactions

% zRiboseInteractions checks all nearby pairs of bases for base-ribose
% interactions, and stores them in a sparse matrix field BaseRibose

function [File,D] = zRiboseInteractions(File,Verbose,Self)

if nargin < 2,
  Verbose = 0;
end

if nargin < 3,
  Self = 1;
end

t = cputime;

% ------------------------------  Meanings of the classification codes:


% ------------------------------ temporary cutoffs! 

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 5.0;                           % near category

AL = 130;                                      % angle limit for BPh
nAL = 90;                                     % angle limit for nBPh

DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;

nDL([4 6 11]) = nNitrogenDist;
nDL([7 8 9])  = nCarbonDist;

% ------------------------------ temporary cutoffs! 

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 5.0;                           % near category

AL = 130;                                      % angle limit for BPh
nAL = 90;                                     % angle limit for nBPh

DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;

nDL([4 6 11]) = nNitrogenDist;
nDL([7 8 9])  = nCarbonDist;

% ------------------------------ specify cutoffs for classification

CarbonDist    = 4.0;                           % max massive - oxygen distance
nCarbonDist   = 4.5;                           % near category

NitrogenDist  = 3.5;                           % max massive - oxygen distance
nNitrogenDist = 4.0;                           % near category

AL = 130;                                      % angle limit for BPh
nAL = 110;                                     % angle limit for nBPh

DL([4 6 11]) = NitrogenDist;
DL([7 8 9])  = CarbonDist;

nDL([4 6 11]) = nNitrogenDist;
nDL([7 8 9])  = nCarbonDist;

% ------------------------------- define basic data

D = [];                          % where to save data if Verbose
T = [];                          % data on which hydrogen with which oxygen(s)

zStandardBases
Sugar = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P','O3 of next'};

p   = [3 5 7];                           % rows of the ribose oxygens
pn  = {'O2*','O3*','O4*'};                % names of ribose oxygens

% ------------------------------ loop through files and classify

for f = 1:length(File),
 if File(f).NumNT > 1,
  if ~isfield(File(f),'Distance'),
    File(f).Distance = [];
  end
  if isempty(File(f).Distance),
    clear c
    for n = 1:length(File(f).NT),
      if ~isempty(File(f).NT(n).Center),
        c(n,:) = File(f).NT(n).Center;
      else
        c(n,:) = [1000*n 0 0];              % pretend, put it far away
      end
    end
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angs
  end

  File(f).BaseRibose = sparse([],[],[],File(f).NumNT,File(f).NumNT);

  % -------- First screening of base pairs ----------------------------------- 

  DistCutoff = 16;                              % max distance for interaction
                                                % screen by C-C distance
  [i,j,v] = find(File(f).Distance);                  % screen by C-C distance
  k = find(v < DistCutoff);
  i = i(k);
  j = j(k);
  v = v(k);



  if Self > 0,
    i = [i; (1:length(File(f).NT))'];             % allow self interactions
    j = [j; (1:length(File(f).NT))'];             % allow self interactions
%    i = [(1:length(File(f).NT))'; i];             % allow self interactions
%    j = [(1:length(File(f).NT))'; i];             % allow self interactions
  end

  % -------- Screen and analyze pairs ----------------------------------------

  for k = 1:length(i),                          % loop through possible pairs
  if ~isempty(File(f).NT(i(k)).Fit) && ~isempty(File(f).NT(j(k)).Fit),

   N1 = File(f).NT(i(k));                       % nucleotide i information
   N2 = File(f).NT(j(k));                       % nucleotide j information

   DT = [];                                     % accumulate data here

   ph = (N2.Sugar(10,:)-N1.Fit(1,:)) * N1.Rot;  % phosphorus displacement
   if abs(ph(3)) < 4.5,                         % phosphorus close to plane

    switch N1.Code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

    dis = zDistance(N1.Fit(m,:), N2.Sugar(p,:)); % distances between mass & O's
    nearDL = nDL(m)' * ones(1,3);     % limits to compare to
    dis = dis .* (dis < nearDL);      % massive-oxygen pairs close enough

    g = [];                           % internal classification number
    w = [];                           % which oxygen interacts
    Angle = [];
    Dist  = [];

    for mm = 1:length(m),             % massive atom to consider
     pp = find(dis(mm,:));            % oxygens close enough to consider
     for n = 1:length(pp),            % loop through potential oxygens
      Angle(n)=zAngle(N1.Fit(m(mm),:),N1.Fit(h(mm),:),N2.Sugar(p(pp(n)),:));
                                      % base massive - hydrogen - oxygen angle
      Dist(n) = dis(mm,pp(n));        % distance
     end

     PAngle = zAngle(N1.Fit(m(mm),:),N1.Fit(h(mm),:),N2.Sugar(10,:));
                                  % base massive - hydrogen - phosphorus angle
     PDist  = zDistance(N1.Fit(m(mm),:),N2.Sugar(10,:));        % distance

     [u,v] = min(-Angle+60*Dist);     % order by quality of potential bond

     for n = 1:length(pp),            % loop through potential oxygens
      if Angle(n) > nAL,              % good enough to be "near" base-phosph

        if ((Angle(n) > AL) && (Dist(n) < DL(m(mm)))) % true BPh pair
          g = [g e(mm)];              % assign a non-near class.
          w = [w pp(n)];              % record which oxygen it is 
          T = [T; [f i(k) j(k) e(mm) pp(n)]];
        else
          g = [g e(mm) + 100];        % > 100 means "near"
          w = [w pp(n)];              % record which oxygen it is
        end

        if Verbose > 1,
          % store information for later display
          ox = (N2.Sugar(p(pp(n)),:)-N1.Fit(1,:)) * N1.Rot; % oxygen displ

          a = [f i(k) j(k) N1.Code g(end) mm pp(n) Angle(n) Dist(n) ox ph File(f).Distance(i(k),j(k)) (v(1)==n) -u(1) PAngle PDist str2num(N1.Number) str2num(N2.Number)];

% [ g(end) v(1) == n]

          % Columns:
          % 1  file number
          % 2  index of base
          % 3  index of nucleotide using ribose
          % 4  code of base
          % 5  classification number for this massive-oxygen pair
          % 6  which massive atom is interacting
          % 7  which oxygen is interacting
          % 8  angle of interaction, in degrees
          % 9  distance from massive atom to oxygen, in Angstroms
          %10  displacement of oxygen atom relative to C1' of base
          %13  displacement of phophorus atom relative to C1' of base
          %16  distance between centers of the two bases
          %17  1 if this is the best oxygen for this hydrogen, 0 otherwise
          %18  approximate quality of the hydrogen bond
          %19  angle made by massive, hydrogen, phosphorus
          %20  distance from massive to phosphorus
          %21  nucleotide number of base
          %22  nucleotide number of ribose donor
          %23  (to be added below) classification of this interaction

          if length(a(1,:)) == 22,
            DT = [DT; a];                  % append data to data matrix
          end

%          zOtherRiboseInteractions

          if Verbose > 3,
            fprintf('%6s base %s%5s %3s BPcode %3d %4s ribose donor %s%5s %13s length %6.2f angle %6.2f interaction %s', File(f).Filename, N1.Base, N1.Number, AtomNames{h(mm),N1.Code}, g(end), zBaseRiboseText(g(end)), N2.Base, N2.Number, Sugar{p(pp(n))}, dis(mm,pp(n)), Angle(n), zEdgeText(File(f).Edge(i(k),j(k))));


            if a(17) == 1,
              fprintf('Best\n');
            else
              fprintf('\n');
            end
          end
        end

      end

     end  % loop over potential oxygens
    end   % loop over massive atoms

     if length(g) > 0,
       if (min(g) < 100) && (max(g) > 100),
         w = w(find(g < 100));
         g = g(find(g < 100));               % remove near classifications
       end
       if length(g) > 1,                     % multiple bonds
         g = sort(g);
         if (g(1) == 6) && (g(end) == 8) && (min(w) < max(w)),
           g = 7;                            % two different oxygens
         elseif (g(1) == 6) && (g(end) == 8) && (min(w) == max(w)),
           g = 18;
         elseif (g(1) == 11) && (g(end) == 13) && (min(w) < max(w)),
           g = 12;                           % two bonds formed
         elseif (g(1) == 11) && (g(end) == 13) && (min(w) == max(w)),
           g = 19;                           % two bonds formed
         elseif (min(g) == max(g)),
           g = g(1);
         elseif (g(end) < 100) && (Verbose > 0),
           fprintf('zRiboseInteractions: Another case to consider: ');
           fprintf('File %s Base index %4d Ribose index %4d Classes ', File(f).Filename, i(k), j(k));
           for m = 1:length(g),
             fprintf('%d ', g(m));
           end
           fprintf('\n');
         end
       end

       File(f).BaseRibose(i(k),j(k)) =   g(1);  % record classification

       if Verbose > 1,
         if ~isempty(DT),
           D = [D; [DT g(1)*ones(length(DT(:,1)),1)]];
         end 
       end

     end

   end    % if vertical displacement of ribose is less than 6 Angstroms
  end     % if Fit is not empty
  end     % loop over nucleotide pairs
 else
  File.BaseRibose = [];
 end      % if File.NumNT > 1
end       % loop over files


if Verbose > 1,
  fprintf('Classifying base-ribose interactions took %8.2f minutes\n', (cputime-t)/60);
  zPhosDisplay(D,0);                     % display parameters graphically
end



if Verbose > 2,
  for i = 1:(length(T(:,1))-1),
    if (T(i,2) == T(i+1,2)) && (T(i,3) == T(i+1,3)),
      if min(T(i,4),T(i+1,4)) == 6 && max(T(i,4),T(i+1,4)) == 8,
        T(i,4)   = 7;
        T(i+1,4) = 7;
        T(i,5)   = 5+4*min(T(i,5),T(i+1,5))+max(T(i,5),T(i+1,5));
        T(i+1,5) = 5+4*min(T(i,5),T(i+1,5))+max(T(i,5),T(i+1,5));
      elseif min(T(i,4),T(i+1,4)) == 11 && max(T(i,4),T(i+1,4)) == 13,
        T(i,4)   = 12;
        T(i+1,4) = 12;
        T(i,5)   = 5+4*min(T(i,5),T(i+1,5))+max(T(i,5),T(i+1,5));
        T(i+1,5) = 5+4*min(T(i,5),T(i+1,5))+max(T(i,5),T(i+1,5));
      end
    end
  end

  [Table,Chi1,P,Labels] = crosstab(T(:,4),T(:,5));
end

return

% File = zAddNTData({'1s72','1j5e','2avy','2aw4','2j01'});
% zRiboseInteractions(File,3);

Res = [];
for f = 1:length(File);
  if isempty(File(f).Info.Resolution),
    Res(f) = 10;
  else
    Res(f) = File(f).Info.Resolution;
  end
end
File = File(find(Res <= 3.0));




i = find(T(:,4) == 12);
Search = [T(i,2) T(i,3) T(i,1)];
xDisplayCandidates(File,Search)
