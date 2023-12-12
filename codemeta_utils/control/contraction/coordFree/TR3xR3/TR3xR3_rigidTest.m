function [ ] = TR3xR3_rigidTest( )
%Test point mass moving in 3D with a reference trajectory
% THIS FUNCTION SHOULD BE CALLED INSIDE A BATCH JOB COMMAND IF USING par_TR3xR3_gainGridSearch_outterLoop
% IE: 'job = batch('TR3xR3_rigidTest','Pool',7)'
% where 7 is the number of instances to run which should correspond to
% number of cores
close all; clear all; clc;
%% Default Parameters
mass = 2; %mass
% x0 = rand(3,1); %inital position
x0 = [1;1;1];
v0 = zeros(3,1); %inital velocity
% B = 1; %beta value for contraction (find using beta_bisection)
e_A = [.6;1]; %minimum eigenvalue of the Hessian (half eVal now due to ref. traj.)
e_B = [.6;1]; %maximum eigenvalue of the Hessian
x1 = zeros(3,1); %final position
v1 = zeros(3,1); %final velocity
TFinal = 10; %final time in s
gradf = @(x,xd) (x-xd)/(sqrt(1+norm(x-xd))^2/2); % define our cost function

%% Test if values will give exponential stability using contraction analysis
% [kd,kv,kp,beta,gridData] = par_TR3xR3_gainGridSearch_outterLoop(e_A, e_B, kd_list, kv_list, kp_list);
[kd,kv,kp,beta] = TR3xR3_fminsearch(e_A,e_B,'showprogress','init_gains',[10;10;10],'maxiter',10);
m_contract = TR3xR3_LMIopt2(e_A, e_B, kv, kd, kp, beta); % Find m's by solving again
% fprintf('max B: %f\n',beta);
% if ( any(any(isnan(m_contract))) )
%     error('Values do not produce a valid contraction region...');
% end

%% Simulate control results (has exponential stability)
control=@(t,x) TR3xR3_rigidControlPD(x, v1, mass, kd, kv, gradf);
closedLoop=@(t,x) TR3xR3_rigidModel(x, x1, mass, control(t,x), gradf, kp);

optsOde=odeset('MaxStep',0.01);

[t,x]=ode45(closedLoop,[0 TFinal],[x0;v0;x0/2], optsOde);

save(['data_TR3xR3_' datestr(datetime,'yyyymmdd_HHMMss') '.mat']);
%% Plot upper bounds using M and beta
% Plot states
figure
plot(t,x,'LineWidth',10)
l1 = legend('x','y','z','$\dot{x}$','$\dot{y}$','$\dot{z}$','$x_r$','$y_r$','$z_r$');
set(l1, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('Position');

% Plot convergence
figure
M = kron(m_contract,eye(3));
for i = 1:length(t)
    x_i = x(i,:)'-[x1;v1;x1];
    d(i) = x_i'*M*x_i;
    d_b(i) = d(1)*exp(-beta*t(i));
end

plot(t,d, 'LineWidth', 10);
hold on
plot(t,d_b, 'LineWidth', 10);
l2=legend('System','Min Convergence: exp(-$\beta t$)');
set(l2, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('Distance to Origin');

end
