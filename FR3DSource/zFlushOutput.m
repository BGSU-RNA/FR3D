% Detect whether running on Octave and if so, flush the standard out so people can see the message

function [void] = zFlushOutput

  if exist('OCTAVE_VERSION') ~= 0,
    fflush(stdout);
  else
  	drawnow;
  end
