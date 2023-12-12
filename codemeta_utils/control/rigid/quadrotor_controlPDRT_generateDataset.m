function quadrotor_controlPDRT_generateDataset(nameSequence,GReference_R,GReference_T,GReference_t)
resetRands();
J=diag([5;5;1]);
m=2;
odeMaxStep=0.001;
%odeMaxStep=0.05;

gainRotationError=80;
gainRotationVelocityError=40;
gainTranslationError=3;
gainTranslationVelocityError=3;

GReference=@(t) rotdrd_interpolationProjected(GReference_t,RT2G(GReference_R,GReference_T),t,'spline');
RReference=@(t) G2R(GReference(t));
TReference=@(t) G2T(GReference(t));
TFinal=max(GReference_t);

optsModel={'inertiaMatrix',J,'mass',m};
optsGains={'gainRotationError',gainRotationError,'gainRotationVelocityError',gainRotationVelocityError,...
    'gainTranslationError',gainTranslationError,'gainTranslationVelocityError',gainTranslationVelocityError};
optsControl=[optsModel optsGains];

w0=0.01*randn(3,1);
v0=0.01*randn(3,1);
x0=rigidDyn_statePackRT(RReference(0),w0,TReference(0),v0);

disturbance=@(x) rigidDyn_disturbance_gravity(x,'mass',m);
control=@(t,x) quadrotor_controlPDRT(x,optsControl{:},'reference',RReference(t),TReference(t),...
    'disturbance',disturbance(x));
closedLoop=@(t,x) quadrotor_model(x,optsModel{:},'input',control(t,x),'disturbance',disturbance(x));

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',odeMaxStep);
[t,x]=rigidDyn_odeEuler(closedLoop,[0 TFinal],x0,optsOde);

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
dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference(t),...
    'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError,...
    'disturbance',disturbance(x));
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

fileName=[nameSequence '.mat'];
fprintf('Data saved to %s\n',fileName)
save(fileName)
