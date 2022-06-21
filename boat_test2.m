clear;clc;close all;
global umax river_length xgoal ygoal
%% Trajectory Optimization and Control of Flying Robot Using Nonlinear MPC
% This code finds the optimal trajectory of a boat crossing a river based
% on the cost function specified in BoatCostFcn2. It also plots the
% resulting throttle and steering inputs, and velocity throughout the
% trajectory. 

%% Model
% The states are given as follows
%
% * x(1) - x inertial coordinate of center of mass
% * x(2) - y inertial coordinate of center of mass
% * x(3) - theta boat heading with respect to the y axis
% * x(4) - xdot velocity of x
% * x(5) - ydot velocity of y
% * x(6) - thetadot angular velocity of theta

%Set some of the constants and parameters for the simulations
Ts = 0.4; %Time step
p = 30; %Prediction control horizon
nx = 6; %number of states
nu = 2; %Number of inputs
river_length = 10;
xgoal = river_length;
ygoal = 10;
umax = -3; %Max current speed, or strength of vortex

% Specify goal state
goal_state = [xgoal;ygoal;0;0;0;0];

% Specify the initial conditions for the boat.
x0 = [0;0;-pi/2;0;0;0];  % boat starts at [0,0], facing north
u0 = zeros(nu,1);        % initial inputs are 0

%% Trajectory optimization
%Create an nlmpc object to find the optimal trajectory and control inputs
ny = 6;
nlobj_tracking = nlmpc(nx,ny,nu);

%% 
% State function
nlobj_tracking.Model.StateFcn = "boat_dynamics";

%%
% Specify the prediction and control horizons and place some constraints on
% the states (within the bounds of the river)
nlobj_tracking.Ts = Ts;
nlobj_tracking.PredictionHorizon = 15;
nlobj_tracking.ControlHorizon = 7;
nlobj_tracking.States(1).Max = river_length + .25;
nlobj_tracking.States(1).Min = -.01;
nlobj_tracking.Optimization.CustomIneqConFcn = "BoatIneqConFunction";
%% 
% Place some weights on the outputs and specify a custom cost function
nlobj_tracking.Weights.ManipulatedVariablesRate = 0.2*ones(1,nu);
nlobj_tracking.Weights.OutputVariables = 5*ones(1,nx);
nlobj_tracking.Optimization.CustomCostFcn = "BoatCostFcn2";

%% 
% Set bounds for the inputs
nlobj_tracking.MV(1).Min = -100;
nlobj_tracking.MV(1).Max = 100;
nlobj_tracking.MV(2).Min = -pi/2;
nlobj_tracking.MV(2).Max = pi/2;

%%Check to see if nlmpc object is valid
validateFcns(nlobj_tracking,x0,u0);

%% Simulation
% Simulate the system for up to 80 steps with specified initial conditions.
Tsteps = 80;        
xHistory = x0';
uHistory = [];
lastMV = zeros(nu,1);

%%Generate figure of river with flow field
npoints = 20;
[xv,yv] = meshgrid(linspace(0,river_length,npoints),linspace(-5,20,npoints));

%Comment or uncomment depending on which flow you choose
[q1,q2] = parabolic_flow(xv,yv);
%[q1,q2] = vortex_flow(xv,yv);

figure(1)
hq = quiver(xv,yv,q1,q2);
xlim([-2 river_length+2]);
title('Boat Trajectory')
hold on
yref = [xgoal,ygoal];
plot(xgoal,ygoal,'go','MarkerSize',20)
plot([0,0],[-5,20],'k','LineWidth',3)
plot([river_length,river_length],[-5,20],'k','LineWidth',3)
options = nlmpcmoveopt;

%Track time and throttle cost totals over trajectory
current_time = 0;
J_tot = 0;

% Use nlmpc move command to find next point
options = nlmpcmoveopt;
for k = 1:Tsteps
    %Current state
    xk = xHistory(k,:)';
    % Compute the optimal control move.
    [uk,options] = nlmpcmove(nlobj_tracking,xk,lastMV,goal_state',[],options);
    % Store the control move and update the last MV for the next step.
    uHistory(k,:) = uk';
    lastMV = uk;
    % Update the real plant states for the next step by solving the
    % continuous-time ODEs based on current states xk and input uk.
    ODEFUN = @(t,xk) boat_dynamics(xk,uk);
    [TOUT,YOUT] = ode45(ODEFUN,[0 Ts], xHistory(k,:)');
    %Plot the trajectory and keep a running total of the cost and time
    figure(1)
    plot(YOUT(:,1),YOUT(:,2),'r')
    xHistory(k+1,:) = YOUT(end,:); 
    current_time = current_time + TOUT(end);
    J_tot = J_tot + (abs(uk(1))*TOUT(end));
    %If the position exceeds the river bounds stop the program
    if xHistory(end,1) > river_length -.25
        fprintf("River Crossed\n")
        break
    end           

end

fprintf('Total Fuel Cost: %.0f\n',J_tot)
fprintf('Time to reach other side: %.2f\n', current_time)

xs = xHistory(:,1);
ys = xHistory(:,2);
error = sqrt(((ys(end)-ygoal)^2)+((xs(end)-xgoal)^2));
fprintf('Error; %.2f\n',error)
T_hist = uHistory(:,1);
times = linspace(0,current_time,length(T_hist));
delta_hist = uHistory(:,2);

%Plot throttle input, steering input, and vessel velocity over trajectory
figure(1)
title('Boat Trajectory')
dim = [.2 .5 .3 .3];
formatSpec = 'River Current = %d.';
str = sprintf(formatSpec,umax);
annotation('textbox',dim,'String',str,'FitBoxToText','on','Position',[2,-7,6,4]);


true_xvel = xHistory(:,4);
true_yvel = xHistory(:,5);

true_vel = sqrt((true_xvel.^2)+(true_yvel.^2));
true_vel = true_vel(1:end-1);

figure(2)
subplot(3,1,1)
plot(times,-uHistory(:,1))
grid on
xlabel('Time (s)')
ylabel('Throttle (m/s/s)')
title('Throttle History')
hold on
subplot(3,1,2)
plot(times,delta_hist)
grid on
xlabel('Time (s)')
ylabel('Steering Angle (rad)')
title('Steering History')
hold on
subplot(3,1,3)
plot(times,true_vel)
grid on
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity History')
hold on


function [q1,q2] = parabolic_flow(xv,yv)
global umax river_length
q1 = 0*xv;
q2 = umax*(1-((xv-(river_length/2))/(river_length/2)).^2);
end

function [q1,q2] = vortex_flow(xv,yv)
global umax
gamma = umax;
Vx = zeros(length(xv),length(yv));
Vy = zeros(length(xv),length(yv));
for i = 1:length(xv)
    for j = 1:length(yv)
        xcoord = xv(i,j);
        ycoord = yv(i,j);
        dx = xcoord - 5;
        dy = ycoord - 5;
        r = sqrt((dx^2)+(dy^2));
        Vx(i,j) = (gamma*dy)/(2*pi*r^2);
        Vy(i,j) = (-gamma*dx)/(2*pi*r^2);
    end
end
q1 = Vx;
q2 = Vy;
end

