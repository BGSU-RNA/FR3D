
function [void] = zLinkFigures(FigList)

      for j=1:length(FigList),
        figure(FigList(j))
        sh(j) = subplot(1,1,1);
        try
	        rotate3d on
	      end
      end

      try
	      linkobj = linkprop(sh,...
	                         {'cameraposition',...
	                          'cameraupvector',...
	                          'cameratarget',...
	                          'cameraviewangle'});
	      set(gcf, 'UserData', linkobj);
	    end