%mSynList.m    %By Ali Mokdad     %April 21, 2006
%This function finds all SYN occurances in a PDB file that has been read.
%It follows the latest definition for SYN (functional SYN) found in my research (see my MS thesis):
%A nucleotide is SYN if the angle CHI between the sugar and the base is such that 0<CHI<+90 degrees
%(instead of the conventional definition for SYN which is -90<CHI<+90 degrees)
%Pay attention that CHI can take values between -180 and +180

%To run this function, do the following
% Filenames='rr0033_5S';
% [File,SIndex]=zAddNTData(Filenames,0);
% [synlist chi_degree]= mSynList(File);

function [synlist, chi_degree]= mSynList(File)

%Declare the main variables:
synlist(1,length(File.NT))=0;
chi_degree(1,length(File.NT))=0;

for i=1:length(File.NT)
  if File.NT(i).Code <= 4,
    %First define the vectors: O4P_C1P; GN_C1P; C1P_GN; C4orC2_GN;
    O4P_C1P     = File.NT(i).Sugar(1,:) - File.NT(i).Sugar(7,:);
    GN_C1P      = File.NT(i).Sugar(1,:) - File.NT(i).Loc(1,:);
    C1P_GN      = - GN_C1P;
    C4orC2_GN   = File.NT(i).Loc(1,:)   - File.NT(i).Loc(2,:);
    
    %%%%%chi angle definition: O4*_C1*_N1_C2 (for pyrimidines) or O4*_C1*_N9_C4 (for purines):
    perp_to_sugar       = cross(O4P_C1P,GN_C1P);
    norm_perp_to_sugar  = norm(perp_to_sugar); 
    perp_to_sugar       = perp_to_sugar/norm_perp_to_sugar;
    
    perp_to_base        = cross(C1P_GN,C4orC2_GN);
    norm_perp_to_base   = norm(perp_to_base);
    perp_to_base        = perp_to_base/norm_perp_to_base;
    
    cross_cross_chi = cross(perp_to_base,perp_to_sugar);
    
    % Take the dot product of the 2 vectors 'perp_to_base' & 'perp_to_sugar' to get cos(chi). Take norm(cross product) to get sin(chi). 
    cos_chi = dot(perp_to_sugar,perp_to_base);
    if dot(cross_cross_chi,C1P_GN) < 0
        sin_chi = norm(cross_cross_chi);
    else
        sin_chi = -norm(cross_cross_chi);
    end
    chi_degree(i) = -180*atan2(sin_chi,cos_chi)/pi; %glycosidic bond angle
    
    % Giving nomenclature according to chi values: anti (most common), or syn
    if (chi_degree(i) >= 0) & (chi_degree(i) <= 90)%This is the new definition for functional syn
        synlist(i) = 1; %SYN
    else
        synlist(i) = 0; %anti
    end
  else
    synlist(i) = NaN;
    chi_degree(i) = NaN;
  end
end
