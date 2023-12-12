function [ ] = contraction3D_rigidTest( )
%Test point mass moving in 3D
close all; clear all; clc;
%% Default Parameters
mass = 2; %mass
% x0 = rand(3,1); %inital position
x0 = [0.12;0.81;0.90];
v0 = zeros(3,1); %inital velocity
kd = 1; %position gain
kv = 2; % velocity gain
% B = 1; %beta value for contraction (find using beta_bisection)
eMin = 0.5; %minimum eigenvalue of the Hessian
eMax = 1; %maximum eigenvalue of the Hessian
x1 = ones(3,1); %final position
v1 = zeros(3,1); %final velocity
TFinal = 10; %final time in s

%% Test if values will give exponential stability using contraction analysis
[m_contract, B] = contraction3D_betaBisection(eMin, eMax, kv, kd);
fprintf('max B: %f\n', B);
[m_contract_LMI, B_LMI] = contraction3D_betaBisection(eMin, eMax, kv, kd,'lmi');
fprintf('max B(LMI): %f\n',B_LMI);
contraction3D_plotEval(eMin,eMax,kv,kd,B_LMI,m_contract_LMI);
if ( isempty(m_contract) )
    error('kv(%f), kd(%f), eMin(%f), eMax(%f) values do not produce a valid contraction region', ...
        kv, kd, eMin, eMax);
end

%% Simulate control results (has exponential stability)
control=@(t,x) contraction3D_rigidControlPD(x, x1, v1, mass, kd, kv);
closedLoop=@(t,x) contraction3D_rigidModel(x, mass, control(t,x));

optsOde=odeset('MaxStep',0.01);

[t,x]=ode45(closedLoop,[0 TFinal],[x0;v0], optsOde);
figure
plot(t,x,'LineWidth',10)
l1 = legend('x','y','z','$\dot{x}$','$\dot{y}$','$\dot{z}$');
set(l1, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('Position');

%% Plot upper bounds using M and beta
figure
M = [eye(3) m_contract(1)*eye(3);m_contract(1)*eye(3) m_contract(2)*eye(3)];
for i = 1:length(t)
    x_i = x(i,:)'-[x1;v1];
    d(i) = x_i'*M*x_i;
    d_b(i) = d(1)*exp(-B*t(i));
end

plot(t,d, 'LineWidth', 10);
hold on
plot(t,d_b, 'LineWidth', 10);
l2=legend('System','Min Convergence: exp(-$\beta t$)');
set(l2, 'Interpreter', 'latex')
xlabel('Time (s)');
ylabel('Distance to Origin');

end
