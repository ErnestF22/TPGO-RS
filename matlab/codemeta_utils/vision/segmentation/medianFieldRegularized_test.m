function medianFieldRegularized_test
resetRands()
lambdaReg=1;
lambda=0.3;
maxIt=100;
sigmaNoise=0.1;

uTruth=repmat([ones(1,15) zeros(1,15)],1,4);
%uTruth=sin(linspace(0,4*pi,100));
Nu=length(uTruth);
k=8;
kHalf=k/2;
idxList=(1:Nu)'*ones(1,2*kHalf)+ones(Nu,1)*[-kHalf:-1 1:kHalf];
idxList=mod(idxList-1,Nu)+1;
kList=k*ones(Nu,1);
disp(idxList)

u=uTruth+sigmaNoise*randn(size(uTruth));
uMedian=medianFieldRegularized(u,idxList,kList,[lambdaReg lambda],...
    'collect','maxIt',maxIt,'rho',0.5);
t=1:Nu;
figure(1)
plot(t,u,'r',t,uMedian(:,end),'b',t,uTruth,'k')
figure(2)
plot(uMedian')

f=zeros(1,maxIt);
for it=1:maxIt
    f(it)=medianFieldRegularizedEnergy(uMedian(:,it),u,idxList,kList,[lambdaReg lambda]);
end
figure(3)
plot(f)

disp('Final energy value')
disp(f(end))
disp('Final energy decrease')
disp(f(end-1)-f(end))
