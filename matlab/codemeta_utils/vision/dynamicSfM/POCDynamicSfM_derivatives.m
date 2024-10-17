function POCDynamicSfM_derivatives
resetRands()
%load file produced by quadrotor_controlPDRT_test_movingReference.m
load('sampleTrajectory')

X=rigidTransform(eye(3),[0;-5;0],randn(3,10));
Xs=X(:,1);
% figure(1)
% plotPoints(T,'-');
% hold on
% plotPoints(X,'x')
% hold off
% axis equal
% xlabel('x');
% ylabel('y');
% zlabel('z');

XCamera=rigidTransform(Rsb,Tsb,X,'references','wc');
Xb=squeeze(XCamera(:,1,:));
% figure(2)
% plotPoints(Xb,'-')
% axis equal

Xbp=multiprodMatVec(Rbs,Xs)-taub;
disp(norm(Xb-Xbp,'fro'))
%dXb=multiprodMatVec(dRbs,Xs)-dtaub;
%dXb=-multiprodMatVec(hwb,multiprodMatVec(Rbs,Xs))+multiprodMatVec(hwb,taub)-nub;
%dXb=-multiprodMatVec(hwb,multiprodMatVec(Rbs,Xs)-taub)-nub;
dXb=-multiprodMatVec(hwb,Xb)-nub;
%funCheckDerInterpInterp(t,Xb,dXb)

hdwb=hat3(dwb);
hwbSq=multiprod(hwb,hwb);
%ddXb=-multiprodMatVec(hdwb,Xb)-multiprodMatVec(hwb,dXb)-dnub;
%ddXb=-multiprodMatVec(hdwb,Xb)+multiprodMatVec(hwbSq,multiprodMatVec(hwb,Xb)+nub)+multiprodMatVec(hwb,nub)-alphab;
ddXb=multiprodMatVec(hwbSq-hdwb,Xb)+2*multiprodMatVec(hwb,nub)-alphab;
funCheckDerInterpInterp(t,dXb,ddXb)



