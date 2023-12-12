function rotDyn_controlPD_test
resetRands();
TFinal=10;
J=diag([5;2;1]);
gainVel=3;

RReference=rot_randn();
optsModel={'inertiaMatrix',J};
optsReference={'RReference',RReference};
optsCost=[optsModel optsReference {'gainRotationError',5}];
optsControl=[optsCost {'flagCancelDynamics',true,'flagVelocityInertiaWeight',false,'gainVelocityError',gainVel}];

R0=diag([-1,-1,1]);
w0=0.02*[1;1;1];

x0=rotDyn_statePack(R0(:),w0);

control=@(t,x) rotDyn_controlPD(x,optsControl{:});
closedLoop=@(t,x) rotDyn_model(x,optsModel{:},'torque',control(t,x));

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
switch 2
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end

x=x';
[R,w]=rotDyn_stateUnpack(x);

figure(2)
subplot(3,1,1)
plot(t,rot2euler(R))
subplot(3,1,2)
[c,gradc]=cost(x,optsCost{:},'methodGradient','packed');
semilogy(t,c)
subplot(3,1,3)
switch 2
    case 1
        [~,nw]=cnormalize(x(10:12,:));
        plot(t,nw)
    case 2
        dx=closedLoop(t,x);
        dc=rotDyn_metricPacked(gradc,dx);
        dcExpected=-gainVel*sum(w.^2);
        plot(t,dc,'-',t,dcExpected,'x')
    case 3
        plotRotationTrajectory(R)
        hold on
        plot3dframe(zeros(3,1),RReference,'ko')
        hold off
end

function [phi,gradphi]=cost(x,varargin)
[R,w]=rotDyn_stateUnpack(x);
[phi,gradphi]=rotDyn_controlPD_cost(R,w,varargin{:});
