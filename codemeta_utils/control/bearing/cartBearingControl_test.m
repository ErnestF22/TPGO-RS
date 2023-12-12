function cartBearingControl_test
seed=3;
resetRands(seed)

L=10;
NX=20;
TFinal=7;
dmax=2;

numDataset=1;
%methodIntegration='DT';
methodIntegration='CT';
dt=0.1;

switch numDataset;
    case 1
        X=L*rand(2,NX);
        %XGoal=[L/2;L/2];
        XGoal=L*rand(2,1);

        x0=[L/2;L/2;0];
        ax=[0 L 0 L];
    case 2
        X=L*rand(2,NX);
        XGoal=[L/2;3/2*L];
        x0=[L/2;-L/2;pi/2];
        ax=[0 L -L 2*L];
end
%funs=bearingCostFunctions('angleSq');
funs=bearingCostFunctions('angleL1L2',0.1);

controlArgs={funs,'dmax',dmax,'maxLinearSpeed',2};

bearingCostDisplayOld(X,XGoal,funs,L)
axis(ax)
pause(0.01)

switch lower(methodIntegration)
    case 'ct'
        [t,x]=ode45(@(t,x) model(x,X,XGoal,controlArgs), [0 TFinal], x0);
        x=x';
    case 'dt'
        t=linspace(0,TFinal,round(TFinal/dt)+1);
        Nt=length(t);
        x=zeros(3,Nt);
        x(:,1)=x0;
        for it=2:Nt
            dx=model(x(:,it-1),X,XGoal,controlArgs);
            x(:,it)=x(:,it-1)+(t(it)-t(it-1))*dx;
        end
end

hold on
plot(x0(1),x0(2),'b*')
plot(x(1,:),x(2,:))
hold off

save([mfilename '_data'])

function dx=model(x,X,XGoal,controlArgs)
[YEval,YGoal]=bearingComputeOld(x(1:2),X,XGoal);
u=cartBearingControl(x(3),YEval,YGoal,controlArgs{:});
dx=cartModel(x,u);
