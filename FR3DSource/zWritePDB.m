% zWritePDB(File) writes the nucleotides in File to a PDB file

function [void] = zWritePDB(File,Filename)

fid = fopen(Filename,'w');       % open for writing

a = 1;                                         % atom number

for n = 1:length(File.NT),
  a = zWriteNucleotidePDB(fid,File.NT(n),a);
end

fclose(fid);
