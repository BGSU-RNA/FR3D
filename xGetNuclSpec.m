% xGetNuclSpec parses the text for nucleotide mask, angle weight, location weight

function [OK,LW,AW] = xGetNuclSpec(str)

OK    = [];
AW    = [];
LW    = [];

if length(str) > 0,

  str    = upper(regexprep(str,';| ',',')); % replace delims by commas
  while strfind(str,',,'),
    str = regexprep(str,',,',',');
  end
  str    = regexprep(str,'LW,','LW');       % replace delims by commas
  str    = regexprep(str,'AW,','AW');       % replace delims by commas
  str    = regexprep(str,'LW','l');         % replace LW by l
  str    = regexprep(str,'AW','a');         % replace AW by a
  commas = strfind(str,',');                % find locations of commas
  lim    = [0 commas length(str)+1];        % locations of tokens
  
  for i=1:length(lim)-1                     % loop through tokens
        
    Token = str(lim(i)+1:lim(i+1)-1);       % extract next token
    if Token(1) == '~',                     % opposite of usual sense
      Reverse = 1;
      Token   = Token(2:length(Token));
    else
      Reverse = 0;
    end
  
    if Token(1) == 'l',                     % location weight
      LW = str2num(Token(2:length(Token))); % extract number
    elseif Token(1) == 'a',                 % angle weight
      AW = str2num(Token(2:length(Token))); % extract number
    else                                    % nucleotide being specified
      Nucl = zeros(1,4);                    % start with 0
      for i=1:length(Token),                % go through letters here
        switch Token(i)
          case 'A', Nucl = Nucl + [1 0 0 0];
          case 'C', Nucl = Nucl + [0 1 0 0];
          case 'G', Nucl = Nucl + [0 0 1 0];
          case 'U', Nucl = Nucl + [0 0 0 1];
          case 'M', Nucl = Nucl + [1 1 0 0];
          case 'R', Nucl = Nucl + [1 0 1 0];
          case 'W', Nucl = Nucl + [1 0 0 1];
          case 'S', Nucl = Nucl + [0 1 1 0];
          case 'Y', Nucl = Nucl + [0 1 0 1];
          case 'K', Nucl = Nucl + [0 0 1 1];
          case 'V', Nucl = Nucl + [1 1 1 0];
          case 'H', Nucl = Nucl + [1 1 0 1];
          case 'D', Nucl = Nucl + [1 0 1 1];
          case 'B', Nucl = Nucl + [0 1 1 1];
          case 'N', Nucl = Nucl + [1 1 1 1];
        end
      end
      Nucl = min([1 1 1 1],Nucl);           
      if max(Nucl) > 0,                         % something specified
        if isempty(OK),                         % no mask yet
          if Reverse == 0,                      % 
            OK = Nucl;
          else
            OK = 1 - Nucl; 
          end
        else
          if Reverse == 0,
            OK = max(OK, Nucl);
          else
            OK = min(OK, 1-Nucl);
          end
        end
      end

    end
  end
end

if isempty(OK),
  OK = [1 1 1 1];
end

if isempty(LW),
  LW = 1;
end

if isempty(AW),
  AW = 1;
end
