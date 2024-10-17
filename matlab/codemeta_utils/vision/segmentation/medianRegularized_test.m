function medianRegularized_test
resetRands()
NPoints=7;
uAll=10*(rand(1,NPoints)-0.5);
uReg=-5;%uAll(1);
u=uAll(2:end);
lambda=0.5;
lambdaReg=1;
uMedian=medianRegularized(u,lambda,uReg,lambdaReg);
disp(uMedian)

f=@(uMedian) medianRegularizedEnergy(uMedian,u,lambda,uReg,lambdaReg);


funPlot(f,linspace(-5,5,100));
hold on
ax=axis();
plot([uMedian uMedian],ax(3:4),'k--')
hold off
