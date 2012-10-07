% xListCandidates(Search) prints a candidate list to the screen
% The optional argument NumToOutput limits the list's length
% The optional argument WheretoOutput has this effect:
%   Value 1 : prints a wide listing to the Matlab command window
%   Value 2 : prints a wide listing to an Editbox
%   Value 3 : prints a narrow listing to an Editbox
% The PC compiled version does both 2 and 3.

% It may be run directly from Matlab using the command:
%   xListCandidates(Search);

function [void] = xListCandidates(Search,NumToOutput,WheretoOutput)

File        = Search.File;

Query       = Search.Query;
Candidates  = Search.Candidates;

[s,t]       = size(Candidates);
N           = Query.NumNT;

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

Text{1} = Search.SaveName;
Text{2} = Search.Query.Name;
Text{3} = Search.Query.Description;

Text{t} = '';

if isfield(Search,'AvgDisc'),
  Text{t} = [Text{t} sprintf('  Filename Avg Discrep ')];
  for i=1:N,
    Text{t} = [Text{t} sprintf('%7d ', i)];
  end
elseif Query.Geometric > 0,
  Text{t} = [Text{t} sprintf('  Filename Discrepancy ')];
  for i=1:N,
    Text{t} = [Text{t} sprintf('%7d ', i)];
  end
else
  Text{t} = [Text{t} sprintf('  Filename Number Nucleotides')];
end

c = 'Chains                                             ';
Text{t} = [Text{t} sprintf('%s', c(1:N))];

if WheretoOutput < 3,
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
end   

if N == 2,
  Text{t} = [Text{t} sprintf(' Pair data')];
end
  
% -------------------------------------- list candidates

Config = {'A' , 'S'};

for i=1:min(s,NumToOutput),

  Text{i+t} = '';

  f = double(Candidates(i,N+1));               % file number for this candidate
  Text{i+t} = [Text{i+t} sprintf('%10s', File(f).Filename)];

  Indices = Candidates(i,1:N);                 % indices of nucleotides

  if isfield(Search,'AvgDisc'),
    Text{i+t} = [Text{i+t} sprintf('%12.4f',Search.AvgDisc(i))];
  elseif Query.Geometric > 0,
    Text{i+t} = [Text{i+t} sprintf('%12.4f',Search.Discrepancy(i))];
  else
    Text{i+t} = [Text{i+t} sprintf('%6d',Search.Discrepancy(i))];      % original candidate number
  end

  for j=1:N,
    Text{i+t} = [Text{i+t} sprintf('%3s',File(f).NT(Indices(j)).Base)];    
    Text{i+t} = [Text{i+t} sprintf('%5s',File(f).NT(Indices(j)).Number)];    
  end

  Text{i+t} = [Text{i+t} sprintf(' ')];

  for j=1:N,
    Text{i+t} = [Text{i+t} sprintf('%s',File(f).NT(Indices(j)).Chain)];
  end

  if WheretoOutput < 3,
    for k=1:length(Indices),
      for j=(k+1):length(Indices),
          Text{i+t} = [Text{i+t} sprintf('%6s', zEdgeText(File(f).Edge(Indices(k),Indices(j))))];
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
  end
    
  if N == 2,                        % special treatment for basepairs
    CP(i) = norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                          File(f).NT(Candidates(i,2)).Sugar(1,:));
    Text{i+t} = [Text{i+t} sprintf('   C1*-C1*: %8.4f', CP(i))];
    NT1 = File(f).NT(Candidates(i,1));
    NT2 = File(f).NT(Candidates(i,2));
    Edge  = full(File(f).Edge(Candidates(i,1),Candidates(i,2)));
    Text{i+t} = [Text{i+t} sprintf(' %s ', zEdgeText(Edge))];
    Text{i+t} = [Text{i+t} sprintf('%7.1f ', Edge)];
    SA = {'A', 'S'};
    Text{i+t} = [Text{i+t} sprintf('%c', SA{1+File(f).NT(Candidates(i,1)).Syn})];
    Text{i+t} = [Text{i+t} sprintf('%c', SA{1+File(f).NT(Candidates(i,2)).Syn})];
  end

end

% -------------------------------------- Additional notifications and info
if (Query.Geometric > 0),
  if (Query.RelCutoff > Query.DiscCutoff) && ~isfield(Search,'AvgDisc'),
    fprintf(fidOUT,'Some motifs with discrepancy between %7.4f and %7.4f might not appear above\n\n', Query.DiscCutoff, Query.RelCutoff);
  end
end

if s > NumToOutput,
  fprintf('Only the first %d candidates were listed.\n', NumToOutput);
end

if N == 2,
  figure
  clf
  hist(CP,30)
  fprintf('Average C1''-C1'' distance is: %8.4f\n', mean(CP));
end

% -------------------------------------- Display the listing

if WheretoOutput == 3,
  mEditbox(Text,'List of Candidates',10);
elseif WheretoOutput == 2,
  mEditbox(Text,'Wide list of Candidates',7);
else
  for i=1:length(Text),
    fprintf('%s\n',Text{i});
  end
end
