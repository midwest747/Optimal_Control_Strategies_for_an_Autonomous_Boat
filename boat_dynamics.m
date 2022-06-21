function dydt = boat_dynamics(y, u)
m = 150; %mass in kg
l = 4; %boat is 2*l m long, assume cg is at center
d = 1; %Distance from cg to center of drag
w = 2; %boat width
h = 1; %submerged depth
I = m*((w^2)+(l^2))/12; %Moment of inertia where we approximate as a rectangular prism
rho = 1; %Density of water in kg/m3
Cd = .1; %Estimated drag coefficient

T = u(1);
delta = u(2);
x = y(1);
yy = y(2);
theta = y(3);
xdot = y(4);
yydot = y(5);
thetadot = y(6);

Areay = w*h*((cos(theta))^2)+2*l*h*((sin(theta))^2);
Areax = abs(Areay*(abs(theta)-(pi/2)));
Fdragx = 0.5*rho*((xdot^2)+(yydot^2))*Areax*Cd*sin(theta);
Fdragy = 0.5*rho*((xdot^2)+(yydot^2))*Areay*Cd*cos(theta);

b1 = (1/m)*(T*sin(theta+delta)+Fdragx);
b2 = (1/m)*(T*cos(theta+delta)+Fdragy);
b3 = -T*inv(I)*l*sin(delta)+inv(I)*d*Fdragx*cos(theta)-inv(I)*d*Fdragy*sin(theta);

dydt(1) = xdot;
dydt(2) = yydot;
dydt(3) = thetadot;
dydt(4) = b1;
dydt(5) = b2;
dydt(6) = b3;
dydt = applydisturbance(dydt,y);
dydt = dydt';
end

function [dydt_dist] = applydisturbance(dydtv,y)
global umax
    %umax = 5;
    %Uncomment for parabolic flow field along y
    u = umax*(1-((y(1)-(10/2))/(10/2))^2);
    dydtv(2) = dydtv(2)+u;

% %Uncomment for vortex flow centered at (5,5)
% gamma = umax;
% yr = y(2)-5;
% xr = y(1)-5;
% r = sqrt((yr^2)+(xr^2));
% ux = gamma*(y(2)-5)/(2*pi*r^2);
% uy = -gamma*(y(1)-5)/(2*pi*r^2);
% dydtv(2) = dydtv(2) + uy;
% dydtv(1) = dydtv(1) + ux;
    
%     %Add wind
%         wind_direction = pi/3; %0 degree wind
%         wind_vel = 2; %Wind speed in m/s
%         dydtv(2) = dydtv(2) + wind_vel*sin(wind_direction);
%         dydtv(1) = dydtv(1) + wind_vel*cos(wind_direction);

    dydt_dist = dydtv;
end