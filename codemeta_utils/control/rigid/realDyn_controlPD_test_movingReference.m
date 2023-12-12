function realDyn_controlPD_test_movingReference
resetRands();
m=2;
gainVel=5;

T0=zeros(3,1);
T1=[1;0;0];
T2=[1;1;0];
T3=[-1;-1;1];
TReference_T=[T0 T1 T2 T3];
TReference_t=[0 5 10 15];
TReference=@(t) interp1(TReference_t,TReference_T',t,'linear')';

optsModel={'mass',m};
optsControl=[{'gainTranslationError',15,'gainVelocityError',gainVel} optsModel];

v0=zeros(3,1);

TFinal=max(TReference_t);

x0=realDyn_statePack(T0,v0);

control=@(t,x) realDyn_controlPD(x,optsControl{:},'TReference',TReference(t));
closedLoop=@(t,x) realDyn_model(x,optsModel{:},'force',control(t,x));

figure(1)
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
[t,x]=ode45(closedLoop,[0 TFinal],x0,optsOde);

x=x';
[T,v]=realDyn_stateUnpack(x);

figure(2)
subplot(2,1,1)
plot(t,T)
hold on
plot(t,TReference(t),':')
hold off

subplot(2,1,2)
plotPoints(T,'-')
hold on
plotPoints(TReference(t),':')
hold off
axis equal
axis([-1 1 -1 1 -1 1])
