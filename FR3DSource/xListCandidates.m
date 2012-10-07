% xListCandidates(Search) prints a candidate list to the screen
% The optional argument NumToOutput limits the list's length
% The optional argument WheretoOutput has this effect:
%   Value 1 : prints a wide listing to the Matlab command window
%   Value 2 : prints a wide listing to an Editbox
%   Value 3 : prints a narrow listing to an Editbox
%   Value 5 : returns a wide listing, doesn't print anything
%   Value 6 : returns the text of the narrow listing
%   Value 7 : narrow listing with information on the organism

% The PC compiled version does both 2 and 3.

% It may be run directly from Matlab using the command:
%   xListCandidates(Search);

function [Text] = xListCandidates(Search,NumToOutput,WheretoOutput,Param)

File        = Search.File;
Candidates  = Search.Candidates;

if ~isfield(Search,'Query'),
  Query.Geometric = 0;
  Query.Name = '';
else
  Query = Search.Query;
end

[s,t]       = size(Candidates);
N           = t-1;

if s == 0,
  fprintf('There are no candidates to list\n');
  return
end

if N == 2,
  CP = zeros(1,s);
end

if nargin < 2,
  NumToOutput = Inf;                    % limit on number printed to screen
end

if (nargin < 3),
  if isdeployed,
    WheretoOutput = 2;
    xListCandidates(Search,NumToOutput,3);
  else
    WheretoOutput = 1;                    % make a wide listing
  end
end

% -------------------------------------- print header line

t = 4;

if isfield(Search,'SaveName'),
  Text{1} = Search.SaveName;
else
  Text{1} = '';
end

if isfield(Search.Query,'Name'),
  Text{2} = Search.Query.Name;
else
  Text{2} = '';
end

if isfield(Search.Query,'Description'),
  Text{3} = Search.Query.Description;
else
  Text{3} = '';
end

Text{t} = '';

if isfield(Search,'AvgDisc'),
  Text{t} = [Text{t} sprintf('  Filename Avg Discrep ')];
elseif Query.Geometric > 0,
  Text{t} = [Text{t} sprintf('  Filename Discrepancy ')];
else
  Text{t} = [Text{t} sprintf('  Filename Number Nucl ')];
end

for i=1:N,
  Text{t} = [Text{t} sprintf('%7d ', i)];
end

c = 'Chains                                             ';
Text{t} = [Text{t} sprintf('%s', c(1:N))];

if isfield(Search,'GroupLabel'),
  Text{t} = [Text{t} ' Group    '];
end

if any(WheretoOutput == [1 2 5]),
  for i=1:N,
    for j=(i+1):N,
      Text{t} = [Text{t} sprintf('%6s', [num2str(i) '-' num2str(j)])];
    end
  end
  
  c = 'Configuration                                      ';
  Text{t} = [Text{t} sprintf(' %s', c(1:N))];
  
  for i=1:N,
    for j=(i+1):N,
      Text{t} = [Text{t} sprintf('%6s', [num2str(i) '-' num2str(j)])];
    end
  end

 if isfield(File,'BaseRibose'),
  for i=1:N,
    for j=1:N,
      if j ~= i,
        Text{t} = [Text{t} sprintf('%6s', [num2str(i) '-' num2str(j)])];
      end
    end
  end
 end

  for i=1:N,
    for j=1:N,
      if j ~= i,
        Text{t} = [Text{t} sprintf('%6s', [num2str(i) '-' num2str(j)])];
      end
    end
  end

  for i=1:N,
    for j=(i+1):N,
      Text{t} = [Text{t} sprintf('%6s', [num2str(i) '-' num2str(j)])];
    end
  end
  
end   

if N == 2,
  Text{t} = [Text{t} sprintf(' Pair data')];
end
  
% -------------------------------------- list candidates

Config = {'A' , 'S'};

for i=1:min(s,NumToOutput),

  f = double(Candidates(i,N+1));               % file number for this candidate
  Indices = Candidates(i,1:N);                 % indices of nucleotides

  Text{i+t} = '';
  Text{i+t} = [Text{i+t} sprintf('%10s', File(f).Filename)];

  if isfield(Search,'DisttoCenter'),
    Text{i+t} = [Text{i+t} sprintf('%12.4f',Search.DisttoCenter(i))];
  elseif Query.Geometric > 0,
    Text{i+t} = [Text{i+t} sprintf('%12.4f',Search.Discrepancy(i))];
  else
    Text{i+t} = [Text{i+t} sprintf('%12d',Search.Discrepancy(i))];      % original candidate number
  end

  for j=1:N,
    Text{i+t} = [Text{i+t} sprintf('%3s',File(f).NT(Indices(j)).Base)];    
    Text{i+t} = [Text{i+t} sprintf('%5s',File(f).NT(Indices(j)).Number)];    
  end

  Text{i+t} = [Text{i+t} ' '];

  for j=1:N,
    Text{i+t} = [Text{i+t} sprintf('%s',File(f).NT(Indices(j)).Chain)];
  end

  if isfield(Search,'GroupLabel'),
    GL = Search.GroupLabel{i};
    Text{i+t} = [Text{i+t} ' ' GL(1:10)];
  end

  if any(WheretoOutput == [1 2 5]),
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        C1 = File(f).NT(Indices(k)).Code;
        C2 = File(f).NT(Indices(j)).Code;
        Text{i+t} = [Text{i+t} sprintf('%6s', zEdgeText(File(f).Edge(Indices(k),Indices(j)),1,C1,C2))];
      end
    end
    
    Text{i+t} = [Text{i+t} sprintf(' ')];
    
    for k=1:length(Indices),
      Text{i+t} = [Text{i+t} sprintf('%c', Config{File(f).NT(Indices(k)).Syn+1})];
    end
    
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        Text{i+t} = [Text{i+t} sprintf('%6d', abs(double(Indices(k))-double(Indices(j))))];
      end
    end

    for k=1:length(Indices),
      for j=1:length(Indices),
        if j ~= k,
         Text{i+t} = [Text{i+t} sprintf('%6s', zBasePhosphateText(File(f).BasePhosphate(Indices(k),Indices(j))))];
        end
      end
    end

   if isfield(File,'BaseRibose'),
    for k=1:length(Indices),
      for j=1:length(Indices),
        if j ~= k,
         Text{i+t} = [Text{i+t} sprintf('%6s', zBaseRiboseText(File(f).BaseRibose(Indices(k),Indices(j))))];
        end
      end
    end
   end

    for k=1:length(Indices),
      for j=(k+1):length(Indices),
        bbc = max(File(f).Backbone(Indices(j),Indices(k)),File(f).Backbone(Indices(k),Indices(j)));
        Text{i+t} = [Text{i+t} sprintf('%6s', zBackboneText(bbc))];
      end
    end

  end
    
  if N == 2,                        % special treatment for basepairs

    CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                          File(f).NT(Candidates(i,2)).Sugar(1,:));
    Text{i+t} = [Text{i+t} sprintf('   C1*-C1*: %8.4f', CP(i))];
    NT1 = File(f).NT(Candidates(i,1));
    NT2 = File(f).NT(Candidates(i,2));
    Edge= full(File(f).Edge(Candidates(i,1),Candidates(i,2)));
    Text{i+t} = [Text{i+t} sprintf('%7.1f ', Edge)];
%    BP  = full(File(f).BasePhosphate(Candidates(i,1),Candidates(i,2)));
%    Text{i+t} = [Text{i+t} sprintf(' %4s ', zBasePhosphateText(BP))];
    SA = {'A', 'S'};
    Text{i+t} = [Text{i+t} sprintf('%c', SA{1+File(f).NT(Candidates(i,1)).Syn})];
    Text{i+t} = [Text{i+t} sprintf('%c', SA{1+File(f).NT(Candidates(i,2)).Syn})];
    if isfield(File,'Crossing'),
      ii = Candidates(i,1);
      jj = Candidates(i,2);
      if (File(f).Crossing(ii,jj) == 0) % && abs(File(f).Edge(ii,jj)) < 15,
        r = ' Nested';
      elseif File(f).Crossing(ii,jj) > 0,
        r = sprintf(' Crossing %4d', full(File(f).Crossing(ii,jj)));
      else
        r = '';
      end
      Text{i+t} = [Text{i+t} r];
    end
  end

  if isfield(File,'Nucl') && (WheretoOutput < 4),
    a = {};
    for j = 1:N,
      if ~isempty(File(f).Nucl(Candidates(i,j)).Motif),
        a = [a File(f).Nucl(Candidates(i,j)).Motif(1).Name];
      end
    end
%    u = unique(a);
    u = a;
    for uu = 1:length(u),
      Text{i+t} = [Text{i+t} ' ' u{uu}];
    end
  end

  if WheretoOutput == 7,
    Text{i+t} = [Text{i+t} ' ' File(f).Info.Source ' | ' File(f).Info.Descriptor];
  end
end

% -------------------------------------- Additional notifications and info
if (Query.Geometric > 0),
  if (Query.RelCutoff > Query.DiscCutoff) && ~isfield(Search,'AvgDisc'),
    L = length(Text);
    Text{L+1} = sprintf('Some motifs with discrepancy between %7.4f and %7.4f might not appear above\n\n', Query.DiscCutoff, Query.RelCutoff);
  end
end

if s > NumToOutput,
  L = length(Text);
  Text{L+1} = sprintf('Only the first %d candidates were listed.\n', NumToOutput);
end

if (N == 2) && (WheretoOutput < 4),
  figure(12)
  clf
  hist(CP,30)
  title('Histogram of C1''-C1'' distance for these pairs');
  fprintf('Average C1''-C1'' distance is: %8.4f\n', mean(CP));
end

% -------------------------------------- Display the listing

if WheretoOutput == 3,
  mEditbox(Text,'List of Candidates',10);
elseif WheretoOutput == 2,
  mEditbox(Text,'Wide list of Candidates',7);
elseif any(WheretoOutput == [1 7]),
  for i=1:length(Text),
    fprintf('%s\n',Text{i});
  end
end
