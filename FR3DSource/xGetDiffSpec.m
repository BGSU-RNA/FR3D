% xGetDiffSpec parses the text for nucleotide differences

function [mindiff,maxdiff,sign] = xGetDiffSpec(str)

mindiff = 1;                                % defaults
maxdiff = Inf;

sign    = 0;                                % default sign of difference

if length(str) > 0,

  str    = regexprep(str,';| ',',');        % replace delims by commas
  while strfind(str,',,'),
    str = regexprep(str,',,',',');          % remove double commas
  end
  str    = regexprep(str,'>=','g');         % replace 
  str    = regexprep(str,'<=','l');         % replace 
  str    = regexprep(str,'==','=');         % replace 
  str    = regexprep(str,'<,','<');         % replace 
  str    = regexprep(str,'>,','>');         % replace 
  str    = regexprep(str,'l,','l');         % replace 
  str    = regexprep(str,'g,','g');         % replace 
  str    = regexprep(str,'=,','=');         % replace 
  str    = regexprep(str,'d','');           % replace d by nothing
  str    = regexprep(str,'D','');           % replace D by nothing

  lt = strfind(str,'>');
  for i=length(lt):-1:1,
   if lt(i)+1 <= length(str),
    if any(str(lt(i)+1) == 'gl=<>'),
      str = [str(1:lt(i)) ',' str(lt(i)+1:length(str))];    
    end
   end
  end

  lt = strfind(str,'<');
  for i=length(lt):-1:1,
   if lt(i)+1 <= length(str),
    if any(str(lt(i)+1) == 'gl=<>'),
      str = [str(1:lt(i)) ',' str(lt(i)+1:length(str))];    
    end
   end
  end

  commas = strfind(str,',');                % find locations of commas
  lim    = [0 commas length(str)+1];        % locations of tokens
  
  for i=1:length(lim)-1                     % loop through tokens
    Token = str(lim(i)+1:lim(i+1)-1);       % extract next token
    if length(Token) > 1,                   % first character is a symbol
      n = str2num(Token(2:length(Token)));    % extract number
      if ~isempty(n),
        switch Token(1)
          case '<', maxdiff = min(maxdiff,n-1); % reduce maxdiff
          case 'l', maxdiff = min(maxdiff,n);   % reduce maxdiff
          case '>', mindiff = max(mindiff,n+1); % increase mindiff
          case 'g', mindiff = max(mindiff,n);   % increase mindiff
          case '=', mindiff = n;
                    maxdiff = n;
        end
      end
    elseif length(Token) > 0,
      switch Token(1)
        case '<', sign = -1;
        case '>', sign =  1;
      end
    end

  end
end

