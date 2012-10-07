% zArcPlot(x,y,color,p1,v1) plots a curved arc connecting points on the unit circle.  The x coordinates of the points are stored in x, the y coordinates in y.  The arc is colored with color, and uses line properties p1 and v1.

function [void] = zArcPlot(x,y,color,Thickness)

if any(color == 'bcrgymk'),
  color = find(color == 'bcrgymk');
end

switch color
case 1,
  h = [0 0 1];                          % blue
case 2,
  h = [0 1 1];                          % cyan
case 3,
  h = [1 0 0];                          % red
case 4,
  h = [0 1 0];                          % green
case 5,
  h = [1 1 0];                          % yellow
case 6,
  h = [1 0 1];                          % magenta
otherwise,
  h = [0 0 0];                          % black
end

p = [x(1) y(1)];
q = [x(2) y(2)];

d = p(1)*q(2) - p(2)*q(1);

if abs(d) > 0.001,
  c(1) = (q(2)*p*p' - p(2)*q*q')/d;
  c(2) = (p(1)*q*q' - q(1)*p*p')/d;
  r = norm(p-c);
  s = atan2(p(2)-c(2),p(1)-c(1));
  t = atan2(q(2)-c(2),q(1)-c(1));
  if t - s > pi,
   s = s + 2*pi;
  elseif s - t > pi,
   t = t + 2*pi;
  end
  u = min(s,t);
  v = max(s,t);
  a = [u:((v-u)/40):v v];
  hold on
  plot(c(1)+r*cos(a),c(2)+r*sin(a),'Color',h,'LineWidth',Thickness);
else
  plot(x,y,'Color',h,'LineWidth',Thickness);
end