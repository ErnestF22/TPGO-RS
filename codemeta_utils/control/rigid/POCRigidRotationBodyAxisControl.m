function POCRigidRotationBodyAxisControl
%resetRands()
e3=[0;0;1];
kYaw=0;

%task='firstTest';
%task='first';
%task='firstDer';
task='second';

%generate reference normal vector (er) and its derivatives
[er,der,dder]=POCRigidRotationBodyAxisControl_reference(0.2,1);


%check_der(@(t) controlRotAndDer(er(t),der(t),dder(t)))
%check_der(@(t) controlRotDerAndDder(er(t),der(t),dder(t)))


R0=rot_randn();

switch lower(task)
    case 'firsttest'
        tTest=rand;
        controlRotationFirstOrder_test(R0,er(tTest),der(tTest))
    case 'first'
        phiAngle=@(R,er) sphere_dist(R*e3,er)^2/2;
        NIt=1000;
        dt=0.01;
        R=R0;
        c=zeros(1,NIt);
        wAll=zeros(3,NIt);
        t=0;
        tAll=zeros(1,NIt);
        for it=1:NIt
            w=controlRotationFirstOrder(R,er(t),der(t),kYaw);
            t=t+dt;
            R=R*rot(dt*w);
            c(it)=phiAngle(R,er(t));
            wAll(:,it)=w;
            tAll(it)=t;
        end
        subplot(2,1,1)
        semilogy(tAll,c)
        subplot(2,1,2)
        plot(tAll,wAll')
    case 'firstder'
        [RTest,~,~,~,wTest]=rot_randGeodFun();
        wt=@(t) controlRotationFirstOrder(RTest(t),er(t),der(t),kYaw);
        dwt=@(t) controlRotationFirstOrder_der(RTest(t),wTest,er(t),der(t),dder(t),kYaw,0);

        check_der(wt,dwt)
    case 'second'
        kW=1;
        TFinal=10;
        J=diag([5;2;1]);
        R0=rot_randn();
        w0=0.2*sphere_randn();

        x0=rotDyn_statePack(R0(:),w0);

        dx=@(t,x) closedLoop(x,J,er(t),der(t),dder(t),kW,kYaw,0);
        optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
        figure(1)
        [t,x]=rotDyn_odeEuler(dx,[0 TFinal],x0,optsOde);
        [R,w]=rotDyn_stateUnpack(x');
        figure(2)
        subplot(2,1,1)
        plot(t,logrot(R))
        hold on
        plot(t,w)
        hold off
        subplot(2,1,2)
        c=evalfunVec(@(y) cost(y(1:12),J,er(y(13)),der(y(13)),kW,kYaw),[x t]');
        plot(t,c)
%         subplot(3,1,3)
%         dc=evalfunVec(@(y) dcost(y(1:12),J,er(y(13)),der(y(13)),dder(y(13)),kW,kYaw,0),[x t]');
%         plot(t,dc)
end



%compare the total derivative of ev with the desired one
function controlRotationFirstOrder_test(R,er,der)
e3=[0;0;1];
ev=R'*er;
devDesired=sphere_log(ev,e3);
w=controlRotationFirstOrder(R,er,der,0);
dev=R'*der+hat(ev)*w;

disp([devDesired dev])

%control law for first order system
function w=controlRotationFirstOrder(R,er,der,kYaw)
flagSO3Embedding=true;
e3=[0;0;1];
ev=R'*er;
if flagSO3Embedding
    RMin=householderRotation3Min(ev,3);
    w=logrot(RMin)+Rpihalf(ev)*R'*der+kYaw*ev;
else
    w=-Rpihalf(ev)*(sphere_log(ev,e3)-R'*der)+kYaw*ev;
end

%derivative of the control law for first order system
function dw=controlRotationFirstOrder_der(R,w,er,der,dder,kYaw,dkYaw)
e3=[0;0;1];
ev=R'*er;
dev=hat(ev)*w+R'*der;
RMin=householderRotation3Min(ev,3);
DRMin=rot(pi*e3)'*householderRotation_DiffMat(ev,3);
dRMinVec=DRMin*dev;

LogRMin=logrot(RMin);
DLogRMin=rot3_logDiff(eye(3),LogRMin,'tangentVec');
dLogRMin=DLogRMin*dRMinVec;

dRder=dRpihalf(ev,dev)*R'*der+Rpihalf(ev)*hat(w)'*R'*der+Rpihalf(ev)*R'*dder;

dYaw=dkYaw*ev+kYaw*dev;

dw=dLogRMin+dRder+dYaw;

%control law for second order system
function Gamma=controlRotationSecondOrder(x,J,er,der,dder,kW,kYaw,dkYaw)
[R,w]=rotDyn_stateUnpack(x);
wd=controlRotationFirstOrder(R,er,der,kYaw);
dwd=controlRotationFirstOrder_der(R,w,er,der,dder,kYaw,dkYaw);
Gamma=hat(w)*J*w+J*dwd-kW*(w-wd);

function dx=closedLoop(x,J,er,der,dder,kW,kYaw,dkYaw)
Gamma=controlRotationSecondOrder(x,J,er,der,dder,kW,kYaw,dkYaw);
dx=rotDyn_model(x,J,Gamma);

function c=cost(x,J,er,der,kW,kYaw)
e3=[0;0;1];
[R,w]=rotDyn_stateUnpack(x);
wd=controlRotationFirstOrder(R,er,der,kYaw);
c=(sphere_dist(R*e3,er)^2+kW*(w-wd)'*J*(w-wd))/2;

function dc=dcost(x,J,er,der,dder,kW,kYaw,dkYaw)
[R,w]=rotDyn_stateUnpack(x);
e3=[0;0;1];
ev=R'*er;
dev=R'*der+hat(ev)*w;
wd=controlRotationFirstOrder(R,er,der,kYaw);
dwd=controlRotationFirstOrder_der(R,w,er,der,dder,kYaw,dkYaw);
Gamma=hat(w)*J*w-(w-wd)+dwd;
dx=rotDyn_model(x,J,Gamma);
[dR,dw]=rotDyn_stateUnpack(dx);
dc=sphere_log(ev,e3)'*dev+kW*(w-wd)'*J*(dw-dwd);
% if dc>0.005
%     keyboard
% end


%rotation of pi/2 around v (assumes norm(v)==1)
function R=Rpihalf(v)
R=hat(v)+v*v';   

function dR=dRpihalf(v,dv)
dR=hat(dv)+dv*v'+v*dv';   

% function [R0,dR0v,ddR0v]=controlReferenceRotation(x,dx,ddx)
% e3=[0;0;1];
% 
% v=x+e3;
% vp=v/norm(v);
% [dvp,ddvp]=normDerAndDDer(v,dx,ddx);
% 
% R0=2*(vp*vp')-eye(3);
% dR0v=-2*hat(vp)*dvp;
% ddR0v=-2*hat(vp)*ddvp;
% 
% function [dxp,ddxp]=normDerAndDDer(x,dx,ddx)
% I=eye(size(x,1));
% nx=norm(x);
% xp=x/nx;
% 
% Dxp=(I-(xp*xp'))/nx;
% dxp=Dxp*dx;
% ddxp=-1/nx*(2*dxp*xp'+xp*dxp')*dx+Dxp*ddx;
% 
% function [R0,dR0]=controlRotAndDer(x,dx,ddx)
% [R0,dR0v,ddR0v]=controlReferenceRotation(x,dx,ddx);
% dR0=R0*hat(dR0v);
% 
% function [dR0v,ddR0v]=controlRotDerAndDder(x,dx,ddx)
% [~,dR0v,ddR0v]=controlReferenceRotation(x,dx,ddx);
