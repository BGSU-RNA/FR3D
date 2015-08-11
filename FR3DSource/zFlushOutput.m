% Detect whether running on Octave and if so, flush the standard out so people can see the message

function [void] = zFlushOutput

  if exist('OCTAVE_VERSION', 'builtin') == 5,
    fflush(stdout);
    warning('off');
  else
  	drawnow;
  end
