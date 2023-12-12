function rigidDyn_controlPD_test
resetRands();
TFinal=10;
J=diag([5;2;1]);
m=1;
gainVel=3;
flagGravity=true;

RReference=rot_randn();
TReference=randn(3,1);
optsModel={'inertiaMatrix',J,'mass',m};
optsReference={'reference',RReference,TReference};
optsCost=[optsModel optsReference {'gainRotationError',5,'gainTranslationError',5}];
optsControl=[optsCost {'gainRotationVelocityError',gainVel,'gainTranslationVelocityError',gainVel}];

R0=diag([-1,-1,1]);
T0=[-1;-1;-1];
w0=0.02*[1;1;1];
v0=0.2*[1;1;1];
x0=rigidDyn_statePackRT(R0,w0,T0,v0);

if ~flagGravity
    control=@(t,x) rigidDyn_controlPDRT(x,optsControl{:});
    closedLoop=@(t,x) rigidDyn_model(x,optsModel{:},'input',control(t,x));
else
    disturbance=@(x) rigidDyn_disturbance_gravity(x,'mass',m);
    control=@(t,x) rigidDyn_controlPDRT(x,optsControl{:},'disturbance',disturbance(x));
    closedLoop=@(t,x) rigidDyn_model(x,optsModel{:},'input',control(t,x),'disturbance', disturbance(x));
end

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
switch 2
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rigidDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end

x=x';
[R,w,T,v]=rigidDyn_stateUnpackRT(x);

figure(2)
subplot(3,2,1)
plot(t,rot_logVec(R))
subplot(3,2,2)
plot(t,T)
subplot(3,1,2)
[c,gradc]=cost(x,optsCost{:},'methodGradient','packed');
semilogy(t,c)
subplot(3,1,3)
switch 4
    case 1
        [~,nw]=cnormalize(x(10:12,:));
        plot(t,nw)
    case 2
        dx=closedLoop(t,x);
        dc=rigidDyn_metricPacked(gradc,dx);
        dcExpected=-gainVel*sum(w.^2)-gainVel*sum(v.^2);
        plot(t,dc,'-',t,dcExpected,'x')
    case 3
        plotRotationTrajectory(R)
        hold on
        plot3dframe(zeros(3,1),RReference,'ko')
        hold off
    case 4
        [dw,dv]=rigidDyn_inputUnpack(control(t,x));
        subplot(3,2,5)
        plot(t,dw)
        subplot(3,2,6)
        plot(t,dv)
end

function [phi,gradphi]=cost(x,varargin)
[R,w,T,v]=rigidDyn_stateUnpackRT(x);
[phi,gradphi]=rigidDyn_controlPDRT_cost(R,w,T,v,varargin{:});
