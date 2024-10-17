function dynamicSfM_integrateIMU_test
load('sampleTrajectory')
[RIntegrated,TIntegrated]=POCIntegrationIMUFunction(t,wIMU,alphaIMU,Rbs(:,:,1),Tsb(:,1),dTsb(:,1),'left');
subplot(2,1,1)
plot(t,reshape(Rbs,9,[]) ,'x',t,reshape(RIntegrated,9,[]),'-')
subplot(2,1,2)
plot(t,Tsb,'x',t,TIntegrated','-')
