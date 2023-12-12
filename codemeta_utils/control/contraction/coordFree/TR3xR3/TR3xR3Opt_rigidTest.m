function [ ] = TR3xR3Opt_rigidTest( )
%Test point mass moving in 3D with a reference trajectory
% THIS FUNCTION SHOULD BE CALLED INSIDE A BATCH JOB COMMAND IF USING par_TR3xR3_gainGridSearch_outterLoop
% IE: 'job = batch('TR3xR3_rigidTest','Pool',7)'
% where 7 is the number of instances to run which should correspond to
% number of cores
close all; clear all; clc;
%% Default Parameters
load('data_TR3xR3_20201019_184715.mat');

%% Generate contraction matrices
M_contract = TR3xR3Opt_contractionMat_general(hessf,hessf, m_contract);
M_contract_der_gains = TR3xR3_contractionMat_general_kdot(hessf,hessf, m_contract);

%% Simulate control results (has exponential stability)
control=@(t,x) TR3xR3Opt_rigidControlPD(x, v1, mass, gradf);
closedLoop=@(t,x) TR3xR3Opt_rigidModel(x, x1, mass, control(t,x), gradf, M_contract, M_contract_der_gains, m_contract);

% optsOde=odeset('MaxStep',0.01);
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);

[t,x]=ode45(closedLoop,[0 TFinal],[x0;v0;x0/2;kd;kv;kp], optsOde);

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
