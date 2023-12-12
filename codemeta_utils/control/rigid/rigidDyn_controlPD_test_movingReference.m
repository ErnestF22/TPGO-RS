function rigidDyn_controlPD_test_movingReference
resetRands();
TFinal=10;
J=diag([5;2;1]);
m=1;
gainVel=3;
gainError=5;

optsModel={'inertiaMatrix',J,'mass',m};
optsControl=[optsModel ...
    {'gainRotationError',gainError,'gainTranslationError',gainError,...
    'gainRotationVelocityError',gainVel,'gainTranslationVelocityError',gainVel}];

T0=zeros(3,1);
T1=[1;0;0];
T2=[1;1;0];
T3=[-1;-1;1];
R0=eye(3);
R1=rot([0;0;1]*pi/2);
R2=R1*rot(-[0;1;0]*pi/2);
R3=R2*rot([0;0;1]*pi/2);
GReference_R=cat(3,R0,R1,R2,R3);
GReference_T=[T0 T1 T2 T3];
GReference_t=[0 5 10 15];
GReference=@(t) rotdrd_interpolationLinear(GReference_t,RT2G(GReference_R,GReference_T),t);
RReference=@(t) G2R(GReference(t));
TReference=@(t) G2T(GReference(t));

w0=zeros(3,1);
v0=zeros(3,1);
x0=rigidDyn_statePackRT(R0,w0,T0,v0);

control=@(t,x) rigidDyn_controlPDRT(x,optsControl{:},'Reference',RReference(t),TReference(t));
closedLoop=@(t,x) rigidDyn_model(x,optsModel{:},'input',control(t,x));

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
subplot(2,2,1)
plot(t,rot_logVec(R),'-',t,rot_logVec(RReference(t)),':')
subplot(2,2,2)
plot(t,T,'-',t,TReference(t),':')
subplot(2,2,3)
plotRotationTrajectory(R)
hold on
plotRotationTrajectory(RReference(t),':')
hold off
subplot(2,2,4)
plotPoints(T,'-')
hold on
plotPoints(TReference(t),':')
