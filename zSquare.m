% zBox(X,Y,color) plots a rectangular box with lower left corner X and
% upper right corner Y, using color

function [void] = zBox(X,Y,color)

set(gcf,'Renderer','OpenGL');

a = X(1);
b = X(2);
c = X(3);

x = Y(1);
y = Y(2);
z = Y(3);

plot3([a,x,x,a,a],[b,b,y,y,b],[c,c,c,c,c],color);
