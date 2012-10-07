% zEdgeFromEdgeText converts a text description such as 'cWW' to a number

function [Edge] = zEdgeFromEdgeText(e)

  e = strrep(e,' ','');

  e = lower(e);

  switch e
    case '' , Edge = 0;
    case 'cww', Edge = 1;
    case 'tww', Edge = 2;
    case 'cwh', Edge = 3;
    case 'twh', Edge = 4;
    case 'cws', Edge = 5;
    case 'tws', Edge = 6;
    case 'chh', Edge = 7;
    case 'thh', Edge = 8;
    case 'chs', Edge = 9;
    case 'ths', Edge = 10;
    case 'css', Edge = 11;
    case 'tss', Edge = 12;
    case 'bif', Edge = 13;
    case 'chw', Edge = -3;
    case 'thw', Edge = -4;
    case 'csw', Edge = -5;
    case 'tsw', Edge = -6;
    case 'csh', Edge = -9;
    case 'tsh', Edge = -10;
    otherwise,  
      Edge = 0;
      fprintf('Unknown edge text in zEdgeFromEdgeText\n');
  end

