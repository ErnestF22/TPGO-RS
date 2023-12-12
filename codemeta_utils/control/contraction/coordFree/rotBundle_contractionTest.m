function [] = rotBundle_contractionTest()
% Test convergence from one rotation to another
close all;

% Default Parameters
kd_list = linspace(.1,100,100);
kv_list = linspace(.1,100,100);
% kd_list = kd_list(end-5:end);
% kv_list = kv_list(end-5:end);
maxD = pi-0.1; % max distance for log(R)
magW = 2; % max angular speed

% Find the best kd, kv, beta
[kd, kv, beta, m_contract, gridData] = rotBundle_gainGridSearch(maxD, kd_list, kv_list, magW);

% Plot the max beta results as function of kd,kv
A = reshape(gridData(:,1),length(kd_list),length(kv_list));
B = reshape(gridData(:,2),length(kd_list),length(kv_list));
C = reshape(gridData(:,3),length(kd_list),length(kv_list));
figure
surf(A,B,C);
xlabel('kd');
ylabel('kv');
zlabel('beta');

% Define the simulation paramters using the best kd, kv, beta
resetRands();
gainVel = kv; % velocity error gain
TFinal=10;
J=diag([5;2;1]);
w0=0.02*cnormalize(randn(3,1));
w = cnormalize(w0); % Choose a random direction
R0 = rot_expVec(eye(3),maxD*w);
w0=-w0;
RReference = eye(3);
optsModel={'inertiaMatrix',J};
optsReference={'RReference',RReference};
optsCost=[optsModel optsReference {'gainRotationError',abs(kd)}];
optsControl=[optsCost {'flagCancelDynamics',true,'flagVelocityInertiaWeight',false,'gainVelocityError',gainVel}];
x0=rotDyn_statePack(R0(:),w0);

% Simulate the system using the given gains
control=@(t,x) rotDyn_controlPD(x,optsControl{:});
closedLoop=@(t,x) rotDyn_model(x,optsModel{:},'torque',control(t,x));

% figure(1)
% optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
optsOde=odeset('MaxStep',0.01);
switch 1
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end
x=x';
[R,w]=rotDyn_stateUnpack(x);
% Save the results for later processing
save(['data_contractionOpt_' datestr(datetime,'yyyymmdd_HHMMss') '.mat']);
end

