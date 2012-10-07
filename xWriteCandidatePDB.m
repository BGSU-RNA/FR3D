% xWriteCandidatePDB(File,Search) writes a single PDB file, named after
% Search, in which each candidate is present, translated from the previous
% one by 20 Angstroms

function [void] = xWriteCandidatePDB(Search)

File = Search.File;

N = Search.Query.NumNT;                        % number of nucleotides in each

if isfield(Search.Query,'LocWeight'),
  LW = Search.Query.LocWeight;
else
  LW = ones(N,1);
end

if isfield(Search.Query,'AngleWeight'),
  AW = Search.Query.AngleWeight;
else
  AW = ones(N,1);
end

if length(Search.Candidates) > 0,

fid = fopen([Search.SaveName '-Cand.pdb'],'w');       % open for writing

M = length(Search.Candidates(:,1));            % number of candidates

f     = Search.Candidates(1,N+1);
Model = File(f).NT(Search.Candidates(1,1:N));  % first cand, taken as model


a = 1;                                         % atom number

for c = 1:M,                                   % loop through candidates
 f     = Search.Candidates(c,N+1);             % file number, this candidate
 Cand  = File(f).NT(Search.Candidates(c,1:N)); % current candidate
 [R,Sh] = xSuperimposeCandidates(Model,Cand,LW,AW);

 for i = 1:N,                                  % loop through nucleotides
  NT = Cand(i);                                % current nucleotide
  a = zWriteNucleotidePDB(fid,NT,a,c,R,Sh);
 end
end

fclose(fid);

fprintf('Wrote %s\n', [Search.SaveName '-Cand.pdb']);

end