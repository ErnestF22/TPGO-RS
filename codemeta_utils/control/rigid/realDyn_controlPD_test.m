function realDyn_controlPD_test
resetRands();
TFinal=10;
m=2;
gainError=5;
gainVel=10;
flagGravity=true;

setup=2;
switch setup
    case 1
        TReference=randn(3,1);
        T0=randn(3,1);
        v0=[1;0;0];
    case 2
        TReference=[0;1;1];
        T0=zeros(3,1);
        v0=zeros(3,1);
end

optsModel={'mass',m};
optsReference={'TReference',TReference};
optsCost=[optsModel optsReference {'gainTranslationError',gainError}];
optsControl=[optsCost {'gainVelocityError',gainVel}];


x0=realDyn_statePack(T0(:),v0);

if ~flagGravity
    control=@(t,x) realDyn_controlPD(x,optsControl{:});
    closedLoop=@(t,x) realDyn_model(x,optsModel{:},'force',control(t,x));
else
    disturbance=@(x) realDyn_disturbance_gravity(x,optsModel{:});
    control=@(t,x) realDyn_controlPD(x,optsControl{:},'disturbance',disturbance(x));
    closedLoop=@(t,x) realDyn_model(x,optsModel{:},'force',control(t,x),'disturbance',disturbance(x));
end

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
[t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);

x=x';
[T,v]=realDyn_stateUnpack(x);

figure(2)
subplot(3,1,1)
plot(t,T)
subplot(3,1,2)
[c,gradc]=cost(x,optsCost{:});
semilogy(t,c)
subplot(3,1,3)
switch 3
    case 1
        nv=cnorm(x(4:6,:));
        plot(t,nv)
    case 2
        dx=closedLoop(t,x);
        dc=sum(gradc.*dx);
        dcExpected=-gainVel*sum(v.^2);
        plot(t,dc,'-',t,dcExpected,'x')
    case 3
        F=control(t,x);
        tInterp=linspace(0,max(t),10);
        TInterp=interp1(t,T',tInterp)';
        FInterp=interp1(t,F',tInterp)';
        plotPoints(T)
        hold on
        plotPoints(TReference,'ko')
        plotField(TInterp,FInterp)
        hold off
        axis equal
end

function [phi,gradphi]=cost(x,varargin)
[T,v]=realDyn_stateUnpack(x);
[phi,gradphi]=realDyn_controlPD_cost(T,v,varargin{:});
