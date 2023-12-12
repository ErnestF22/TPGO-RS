function POCRigidRotationControlReference
resetRands();

%J=eye(3);
J=diag([5;2;1]);

R1=eye(3);
R2=rot([0;0;1]*pi/2);
R3=R2*rot(-[0;1;0]*pi/2);
R4=R3*rot([0;0;1]*pi/2);

RReference_t=[0 5 10 15 20];
%RReference_R=cat(3,R1,R2,R2,R3,R3);
RReference_R=cat(3,R1,R2,R3,R4,R1);
RReference=@(t) rot_interpolationLinear(RReference_t,RReference_R,t);

R0=R1;
w0=0.0*[1;1;1];
TFinal=max(RReference_t);


optsModel={'inertiaMatrix',J};
optsControl={'inertiaMatrix',J,'flagCancelDynamics',true,'gainRotationError',5,'gainVelocityError',5};

x0=rotDyn_statePack(R0(:),w0);

control=@(t,x) rotDyn_controlPD(x,optsControl{:},'RReference',RReference(t));
closedLoop=@(t,x) rotDyn_model(x,optsModel{:},'torque',control(t,x));

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
switch 2
    case 1
        [t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=rotDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);
end

[R,w]=rotDyn_stateUnpack(x');

figure(2)
subplot(2,1,1)
plotRotationTrajectory(R)
hold on
plotRotationTrajectory(RReference(t),'x')
hold off
subplot(2,1,2)
plot(rot_dist(R,RReference(t),'vector'))

