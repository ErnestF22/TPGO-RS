% Script to test Double Integrator System on a Plane
% Combine CLF and CBF to control dbl integrator with safety to follow
% trajectory

close all; clear all; clc;

% Parameters (ASSUME field is [-50,50,-50,50]
m = 10; % Mass
TFinal = 40;
global U_global;
global USE_RRT;
USE_RRT = false;
global path_RRT; % List to hold the waypoints
path_RRT = [];
global idx_RRT; % Index of current waypoint
idx_RRT = 0;
global RRT_StartTime;
RRT_StartTime = 0;
global RRT_StopTime;
RRT_StopTime = 0;
Boundary = [-50;50;-50;50]; % Boundary [xmin,xmax,ymin,ymax]

% Setup and run simulation
[x0, xd, xd_dd, LfB, LgB, Bf, h_i, x_i, y_i, r_i, axisBound] = obs_dynamicTraj(m);
Obstacles =@(t) [x_i(t)-r_i(t),y_i(t)-r_i(t),2*r_i(t),2*r_i(t)];
control = @(t,x) controller(x, xd(t), xd_dd(t), m, LfB(t,x), LgB(t,x),...
    Bf(t,x), h_i(t,x), Boundary, Obstacles(t), t);
closedLoop = @(t, x) model(x, control(t,x), m);
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
figure;
[t,x]=ode45(closedLoop,[0 TFinal],x0, optsOde);
% title('States');
lgd = legend('x','y','dx','dy');
lgd.FontSize = 34;
xlabel('Time (s)','FontSize',20);
ylabel({'(m)','(m/s)'},'FontSize',20);

% Plot trajectory
filename = 'Trajectory.gif';
% figure('units','normalized','outerposition',[0 0 1 1]);
figure
timeStepSize = 0.25; % Time step to simulate starting from 0
idx_Time_Current = find(t > RRT_StartTime,1); % Current index corresponding to t
idx_Time_Prev = idx_Time_Current; % Last index corresponding to t-timeStepSize
xd_traj = xd(t')'; % Get the desired trajectory
% for i = 1:gifStepSize:length(x)
bStartNewGif = true;
while (true)   
    cla
    axis(axisBound)
    xlabel('x (m)','FontSize',20);
    ylabel('y (m)','FontSize',20);
    obstacles = Obstacles(t(idx_Time_Current));
    for j = 1:size(obstacles,1)
        obs = obstacles(j,:);
        rectangle('Position',obs,'FaceColor',[0 .5 .5],'Curvature',[1 1])
        hold on
    end
    % Plot RRT* here
    if (t(idx_Time_Current) > RRT_StartTime && ~isempty(path_RRT))
        plot(path_RRT(:,1),path_RRT(:,2),'r--','LineWidth',5);
        plot(path_RRT(:,1),path_RRT(:,2),'rx','LineWidth',5);
        hold on
    end
    % Plot desired trajectory
    plot(xd_traj(idx_Time_Prev:idx_Time_Current,1),...
        xd_traj(idx_Time_Prev:idx_Time_Current,2),'y-','LineWidth',3);
    rectangle('Position',[xd_traj(idx_Time_Current,1)-0.5 xd_traj(idx_Time_Current,2)-0.5 1 1],'FaceColor',.8242*ones(1,3),'Curvature',[1 1]);
    
    % Plot trajectory
    plot(x(1:idx_Time_Current,1),x(1:idx_Time_Current,2),'b--','LineWidth',3);
    % Plot a tail
    plot(x(idx_Time_Prev:idx_Time_Current,1),x(idx_Time_Prev:idx_Time_Current,2),'g-','LineWidth',3);
    % Plot the point mass
    rectangle('Position',[x(idx_Time_Current,1)-0.5 x(idx_Time_Current,2)-0.5 1 1],'FaceColor','k','Curvature',[1 1]);
    drawnow
    frame = getframe();
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if (bStartNewGif)
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
        bStartNewGif = false;
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    % Find the next time index
    idx_Time_Prev = idx_Time_Current;
    idx_Time_Current = find(t>=t(idx_Time_Current)+timeStepSize,1);
    if ( isempty(idx_Time_Current) )
        break;
    end
end
% comet(x(:,1), x(:,2));
% Plot relaxation term u(3)
figure('units','normalized','outerposition',[0 0 1 1]);
plot(U_global(3,:));
title('Relaxtion Term');
% Plot h_i
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(t)
    h_i_results(i,:) = h_i(t(i),x(i,:));
end
plot(t,h_i_results,'LineWidth',5);
xlabel('Time','FontSize',20);
ylabel('h','FontSize',20);
lgd = legend('Obs 1','Obs 2','Obs3');
lgd.FontSize = 34;
save('TestRun.mat'); % Save everything for use later

function [dx] = model(x, u, m)
% INPUTS:
%   x = state vector [x, y, vx, vy]
%   u = control input into controlled vehicle
% OUTPUTS:
%   dx = new states

% System dynamics
A = [zeros(2,2) eye(2);zeros(2, 4)];
B = [zeros(2,2);eye(2)];

dx = A*x + 1/m*B*u(1:2);

% Store u
global U_global;
U_global = [U_global, u];
end

function [u] = controller(x, xd, xd_dd, m, LfB, LgB, Bf, h_i, boundary, obstacles, t)
% INPUTS:
%   x = state vector [x, y, vx, vy]
%   xd = desired trajectory [xd, yd, vxd, vyd]
%   m = mass in kg
%   LfB = evaulated Lie Der. of B wrt f
%   LgB = evaulated Lie Der. of B wrt g
% OUTPUTS:
%   u = controller [u1; u2]
STOP_TOL = 1e-1; % Min norm squared that RRT* should take over
DIST_TOL = 1; % Min dist to say we're closed to desired location
ANGLE_TOL = 45/180*pi; % Angle tolerance to use original trajectory
MIN_DIST_RESUME_TRAJ = 1e-2; % Min dist to travel away from local min
global USE_RRT;
global path_RRT; % List to hold the waypoints
global idx_RRT; % Index of current waypoint
global RRT_StartTime;
global RRT_StopTime;

if ( USE_RRT )
    xd_temp = path_RRT(idx_RRT,:);
    % Find control of the original trajectory
    u_original = optimizeU(x, xd, xd_dd, m, LfB, LgB, h_i);
    % If control to original trajectory is not going towards local minimum,
    % then resume with original trajectory
    vec2LocalMin = path_RRT(1,:)' - x(1:2);
    angle = acos(dot(u_original(1:2), vec2LocalMin)/(norm(u_original(1:2))*norm(vec2LocalMin)));
    if ( norm(vec2LocalMin) > MIN_DIST_RESUME_TRAJ && ...
            angle > ANGLE_TOL && angle < 2*pi-ANGLE_TOL )
        USE_RRT = false;
        RRT_StopTime = t;
        u = u_original;
        fprintf('Using desired trajectory...\n');
    else
        if (norm(xd_temp'-x(1:2)) < DIST_TOL)
            % Reached the waypoint, go to next
            idx_RRT = idx_RRT+1;
            if (idx_RRT > size(path_RRT,1))
                idx_RRT = size(path_RRT,1);
                % Reached stationary trajectory
                USE_RRT = false;
                RRT_StopTime = t;
                fprintf('End of RRT*...\n');
    %         if (idx_RRT > 3)
    %             USE_RRT = false;
            end
            xd_temp = path_RRT(idx_RRT,:);
            fprintf('Setting waypoint to (%0.2f, %0.2f)...\n', xd_temp(1), xd_temp(2));
        end
        xd = [xd_temp';0;0];
        % Find control of temp waypoint
        u = optimizeU(x, xd, xd_dd, m, LfB, LgB, h_i);
    end
else
    u = optimizeU(x, xd, xd_dd, m, LfB, LgB, h_i);
end

% u
if (isempty(u))
    Bf
    u
    h_i
    x
    % Bdot = LfB + LgB*u(1:2) - gamma_B/Bf
    % Vdot = phi0 + phi1*u(1:2)

    error('Problem is not feasible');
end

% If no control input, not moving (no velocity, and have not reached goal,
% do RRT*
if ( (norm(u(1:2)) < STOP_TOL) && (norm(x(3:4)) < STOP_TOL) && (norm(x(1:2) - xd(1:2)) > DIST_TOL))% && ~USE_RRT)
    % Controller is stuck in local minima, use RRT* to path plan to xd
    fprintf('Using RRT* for trajectory...\n');
    USE_RRT = true;
    idx_RRT = 2; % This should be the waypoint after our current position
    path_RRT = func_RRTStar(x(1:2)', xd(1:2)', boundary, ...
        obstacles, 2, 500);
    figure(1); % Assumes this is the main plot (for now)
    xd_temp = path_RRT(idx_RRT,:);
    fprintf('Setting waypoint to (%0.2f, %0.2f)...\n', xd_temp(1), xd_temp(2));
    RRT_StartTime = t;
end

end

function [u] = optimizeU(x, xd, xd_dd, m, LfB, LgB, h_i)
% INPUTS:
%   x = state vector [x, y, vx, vy]
%   xd = desired trajectory [xd, yd, vxd, vyd]
%   m = mass in kg
%   LfB = evaulated Lie Der. of B wrt f
%   LgB = evaulated Lie Der. of B wrt g
% OUTPUTS:
%   u = controller [u1; u2]

% Parameters
K = 6.34; % Gain for ex*ex term of V
eps = 5.4; % Gain for ex*ev term of V
eta = 0.5; % Exponential convergence rate of V
p_sc = 1e3; % relaxation term compensation
gamma_B = 1e-2; % max growth rate of Bdot leq gamma_B
U_Bounds = [500;500;1e10];
% CLF
ex = x(1:2) - xd(1:2);
ev = x(3:4) - xd(3:4);
V = 1/2*m*(ev'*ev) + 1/2*K*(ex'*ex) + eps*(ex'*ev);
% phi0 = (K*ex'+eps*ev')*x(1:2)+eta*V;
% phi1 = ev'+eps/m*ex';
phi0 = K*ex'*ev + eps*(ev'*ev) - (m*ev'+eps*ex')*xd_dd + eta*V;
phi1 = ev'+eps/m*ex';

% CBF (Calculated in main function as function handles)

% Quadratic Parameters
A_clf = [phi1, -1];
b_clf = -phi0;
A_cbf = [LgB, zeros(size(LgB,1),1)];
b_cbf = -LfB + gamma_B*h_i.^3;
% Matrix to use in solver
A = [A_clf; A_cbf];
b = [b_clf; b_cbf];
H_acc = 2*[1/m^2 0 0;0 1/m^2 0;0 0 p_sc];
F_acc = -2*[0;0;0];

% Solve quadratic program using quadprog
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off','MaxIterations',50);
[u,fval,exitflag,output,lambda] = ...
   quadprog(H_acc,F_acc,A, b,...
   [],[],-U_Bounds,U_Bounds,[],options);
end