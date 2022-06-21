function J = BoatCostFcn2(x,u,e,data)
global xgoal ygoal
Rt = 1;
Rs = .1;
Qx = 50;
Qy = 50;
N = length(x(:,1));
Qv = 0;
v = sqrt(x(:,4).^2 + x(:,5).^2);

xdiff = x(:,1)-xgoal;
ydiff = x(:,2)-ygoal;
if xdiff(end) < 2 && ydiff(end) < 2
    Qx = 300;
    Qy = 300;    
    Qv = 100;
end

J = (xdiff'*Qx*xdiff + ydiff'*Qy*ydiff + v'*Qv*v + u(:,1)'*Rt*u(:,1) + u(:,2)'*Rs*u(:,2))/N ;