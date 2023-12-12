function quadrotor_controlPDRT_test_movingReference
resetRands();
J=diag([5;2;1]);
m=2;
flagGravity=true;
odeMaxStep=0.001;
%odeMaxStep=0.05;

if ~flagGravity
    gainRotationError=65;
    gainRotationVelocityError=15;
    gainTranslationError=1.5;
    gainTranslationVelocityError=1.5;
else
    gainRotationError=80;
    gainRotationVelocityError=40;
    gainTranslationError=3;
    gainTranslationVelocityError=3;
end

T0=[0;-1;1];
T1=[1;0;0];
T2=[-1;1;0];
%T2=[1;1;0];
T3=T2;%[1;-1;1];
T4=[1;-1;1];
R0=eye(3);
R1=rot([0;0;1]*10/180*pi);
R2=rot(-[0;0;1]*10/180*pi);
GReference_R=cat(3,R0,R1,R2,R0,R0);
GReference_T=[T0 T1 T2 T3 T4]/2;
GReference_t=[0 5 10 15 20]/2;
GReference=@(t) rotdrd_interpolationProjected(GReference_t,RT2G(GReference_R,GReference_T),t,'spline');
RReference=@(t) G2R(GReference(t));
TReference=@(t) G2T(GReference(t));
%TFinal=max(GReference_t);
TFinal=7.5;

optsModel={'inertiaMatrix',J,'mass',m};
optsGains={'gainRotationError',gainRotationError,'gainRotationVelocityError',gainRotationVelocityError,...
    'gainTranslationError',gainTranslationError,'gainTranslationVelocityError',gainTranslationVelocityError};
optsControl=[optsModel optsGains];

w0=0.01*randn(3,1);
v0=0.01*randn(3,1);
x0=rigidDyn_statePackRT(R0,w0,T0,v0);

if ~flagGravity
    control=@(t,x) quadrotor_controlPDRT(x,optsControl{:},'reference',RReference(t),TReference(t));
    closedLoop=@(t,x) quadrotor_model(x,optsModel{:},'input',control(t,x));
else
    disturbance=@(x) rigidDyn_disturbance_gravity(x,'mass',m);
    control=@(t,x) quadrotor_controlPDRT(x,optsControl{:},'reference',RReference(t),TReference(t),...
        'disturbance',disturbance(x));
    closedLoop=@(t,x) quadrotor_model(x,optsModel{:},'input',control(t,x),'disturbance',disturbance(x));
end    
figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',odeMaxStep);
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
title('R')
subplot(2,2,2)
plot(t,T,'-',t,TReference(t),':');
title('T')
subplot(2,1,2)
if ~flagGravity
    dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference(t),...
        'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError);
else
    dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference(t),...
        'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError,...
        'disturbance',disturbance(x));
end
plotPoints(T,'-')
hold on
plotPoints(TReference(t),':')
tInterp=linspace(0,max(t),20);
TInterp=interp1(t,T',tInterp)';
dvInterp=interp1(t,dv',tInterp)';
e3Interp=interp1(t,squeeze(R(:,3,:))',tInterp)';
plotField(TInterp,dvInterp,'m');
plotField(TInterp,e3Interp,'b');
hold off
axis equal
view(45,45)

disp([x(:,1), closedLoop(0,x(:,1))])
disp(dv(:,1))

save([mfilename '_data'])
