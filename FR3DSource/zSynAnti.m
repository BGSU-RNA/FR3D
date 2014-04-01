% zSynAnti.m by CLZ 2013-11-06

% This program modifies the approach by Ali Mokdad
% mSynList.m    %By Ali Mokdad     %April 21, 2006
% Pay attention that CHI can take values between -180 and +180

function [File, chi_degree] = zSynAnti(File)

chi_degree(1,length(File.NT))=0;

for i=1:length(File.NT)
  if File.NT(i).Code <= 4,
    %First define the vectors: O4P_C1P; GN_C1P; C1P_GN; C4orC2_GN;
    O4P_C1P     = File.NT(i).Sugar(1,:) - File.NT(i).Sugar(7,:);
    GN_C1P      = File.NT(i).Sugar(1,:) - File.NT(i).Fit(1,:);
    C1P_GN      = - GN_C1P;
    C4orC2_GN   = File.NT(i).Fit(1,:)   - File.NT(i).Fit(2,:);
    
    %%%%%chi angle definition: O4*_C1*_N1_C2 (for pyrimidines) or O4*_C1*_N9_C4 (for purines):
    perp_to_sugar       = cross(O4P_C1P,GN_C1P);
    norm_perp_to_sugar  = norm(perp_to_sugar); 
    if norm_perp_to_sugar ~= 0,
        perp_to_sugar   = perp_to_sugar/norm_perp_to_sugar;
    end
    
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
    % sign of chi_degree matches Bevilacqua 2011 paper on syn and anti
    % this definition matches the IUPAC definition from http://www.chem.qmul.ac.uk/iupac/misc/pnuc2.html#230
    chi_degree(i) = 180*atan2(sin_chi,cos_chi)/pi; %glycosidic bond angle
    
    % Giving nomenclature according to chi values: anti (most common), or syn

    File.NT(i).chi = chi_degree(i);

    if chi_degree(i) > -90 && chi_degree(i) < -45,
        classification = 'Intermediate Syn';
    elseif chi_degree(i) >= -45 && chi_degree(i) < 90,
        classification = 'Syn';
    else
        classification = 'Anti';
    end

    File.NT(i).glycosidicbondorientation = classification;

  else
    chi_degree(i) = NaN;
    File.NT(i).chi = NaN;
    File.NT(i).glycosidicbondorientation = '';
  end
end


