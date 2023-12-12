function quadrotor_controlPDRT_test
resetRands();
TFinal=10;
J=diag([5;2;1]);
m=2;
flagGravity=true;

if ~flagGravity
    gainRotationError=65;
    gainRotationVelocityError=15;
    gainTranslationError=1.5;
    gainTranslationVelocityError=3;
else
    gainRotationError=60;
    gainRotationVelocityError=30;
    gainTranslationError=1.5;
    gainTranslationVelocityError=3;
end

RReference=eye(3);
TReference=[1;1;1];
optsModel={'inertiaMatrix',J,'mass',m};
optsReference={'reference',RReference,TReference};
optsGains={'gainRotationError',gainRotationError,'gainRotationVelocityError',gainRotationVelocityError,...
    'gainTranslationError',gainTranslationError,'gainTranslationVelocityError',gainTranslationVelocityError};
optsControl=[optsModel optsReference optsGains];

R0=eye(3);
w0=zeros(3,1);
T0=zeros(3,1);
v0=zeros(3,1);
x0=rigidDyn_statePackRT(R0,w0,T0,v0);

if ~flagGravity
    control=@(t,x) quadrotor_controlPDRT(x,optsControl{:});
    closedLoop=@(t,x) quadrotor_model(x,optsModel{:},'input',control(t,x));
else
    disturbance=@(x) rigidDyn_disturbance_gravity(x,'mass',m);
    control=@(t,x) quadrotor_controlPDRT(x,optsControl{:},'disturbance',disturbance(x));
    closedLoop=@(t,x) quadrotor_model(x,optsModel{:},'input',control(t,x),'disturbance',disturbance(x));
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
subplot(2,2,1)
plot(t,rot_logVec(R))
title('R')
subplot(2,2,2)
plot(t,T)
title('T')
subplot(2,1,2)
if ~flagGravity
    dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference,...
        'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError);
else
    dv=realDyn_controlPD(x, 'mass',m,'TReference',TReference,...
        'gainTranslationError',gainTranslationError,'gainVelocityError',gainTranslationVelocityError,...
        'disturbance',disturbance(x));
end
plotPoints(T,'-')
hold on
plotPoints(TReference,'ko')
tInterp=linspace(0,max(t),10);
TInterp=interp1(t,T',tInterp)';
dvInterp=interp1(t,dv',tInterp)';
e3Interp=interp1(t,squeeze(R(:,3,:))',tInterp)';
plotField(TInterp,dvInterp,'m');
plotField(TInterp,e3Interp,'b');
hold off
axis equal

disp([x(:,1), closedLoop(0,x(:,1))])
disp(dv(:,1))
