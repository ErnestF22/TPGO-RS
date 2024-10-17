function rigidTransform_test
X=randn(3,1);
methodPoses='reference';
%methodPoses='pose';
[Rt,vt,R0,v0,vVec]=rot_randGeodFun(eye(3));
[Tt,dTt,T0,dT]=real_randGeodFun(randn(3,1));
figure(1)
check_der(@(t) evaluateFandDF(X,Rt,Tt,vVec,dT,methodPoses,t),'function')

wFix=randn(6,1);

figure(2)
check_der(@(t) evaluateDFandDdF(X,Rt,Tt,vVec,dT,methodPoses,wFix,t),'function')


function [f,df]=evaluateFandDF(X,Rt,Tt,vVec,dT,methodPoses,t)
[f,Jf]=rigidTransform(Rt(t),Tt(t),X,'methodAbsolutePoses',methodPoses);
df=Jf*[vVec;dT];


function [df,ddf]=evaluateDFandDdF(X,Rt,Tt,vVec,dT,methodPoses,wFix,t)
[~,Jf,Hf]=rigidTransform(Rt(t),Tt(t),X,'methodAbsolutePoses',methodPoses);
w=[vVec;dT];
df=Jf*wFix;
ddf=zeros(3,1);
for k=1:3
    ddf(k)=w'*Hf(:,:,k)*wFix;
end
