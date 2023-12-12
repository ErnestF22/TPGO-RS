function POCRigidRotationControl
resetRands();
TFinal=10;
J=diag([5;2;1]);
gainVel=5;

RReference=eye(3);
%RReference=rot_randn();
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
        dx=rotDyn_controlPD(x,optsControl{:},'flagExtendedVector',true);
        dc=rotDyn_metricPacked(gradc,dx);
        dcExpected=-gainVel*sum(w.^2);
        plot(t,dc,'-',t,dcExpected,'x')
    case 3
        plotRotationTrajectory(R)
        hold on
        plot3dframe(zeros(3,1),RReference,'ko')
end

function [phi,gradphi]=cost(x,varargin)
[R,w]=rotDyn_stateUnpack(x);
[phi,gradphi]=rotDyn_controlPD_cost(R,w,varargin{:});


% Nt=length(t);
% for it=1:Nt
%     T=zeros(3,1);
%     draw3dcameraFromRT(R(:,:,it),T)
%     axis equal
%     axis([-2 2 -2 2 -2 2]);
%     getframe();
%     pause(0.01)
% end


function dPhi=derPhi(x,J)
gradPhi=gradCost(x,J);
[R,wx]=rotDyn_stateUnpack(x);
dx=closedLoop(x,J);
[dR,dw]=rotDyn_stateUnpack(dx);
w=rot_vee(R,dR);
dPhi=[w;dw]'*gradPhi;


function dx=closedLoop(t,x,varargin)
control=@(x,J) zeros(3,1);
optsModel={};

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'control'
            ivarargin=ivarargin+1;
            control=varargin{ivarargin};
        case 'optsmodel'
            ivarargin=ivarargin+1;
            optsModel=[optsModel varargin{ivarargin}];
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

dx=rotDyn_model(x,optsModel{:},'torque',control(t,x));

% function c=cost(x,J,varargin)
% [R,w]=rotDyn_stateUnpack(x);
% RReference=eye(3);
% 
% %optional parameters
% ivarargin=1;
% while ivarargin<=length(varargin)
%     switch lower(varargin{ivarargin})
%         case 'rreference'
%             ivarargin=ivarargin+1;
%             RReference=varargin{ivarargin};
%         otherwise
%             error('Argument not valid!')
%     end
%     ivarargin=ivarargin+1;
% end
% 
% c=(rot_dist(R,RReference)^2+w'*J*w)/2;
% 
% % function c=gradCost(x,J)
% % [R,w]=rotDyn_stateUnpack(x);
% % 
% % c=[logrot(R);J*w];

function gamma=control(x,J)
[R,w]=rotDyn_stateUnpack(x);
skew=@(A) (A-A')/2;
eR=logrot(R);
%eR=eR/norm(eR);
delta=0.1;
%uw=ew;
ew=w;
[ewNorm,new]=cnormalize(ew);
k1=0.01^2;
if new^2<k1
    uw=ew;
else
    uw=k1*ew/new^2;
end
uw=1000*uw;
%uw=-0.02*cnormalize(eR);
ew=ew+hat3(w)*J*w;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
%ewNorm=norm(ew);
%if ewNorm>delta
%    ew=ew/ewNorm;
%else
%    ew=ew/delta;
%end
gamma=(-eR-uw);

function gamma=control2(x,J)
[R,w]=rotDyn_stateUnpack(x);
eR=logrot(R);
wr=cnormalize(logrot(R));
gamma=hat3(w)*J*w+w-wr;
