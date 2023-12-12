function bearingCostGeneral_desiredAngle_derBound_test
global funs;


%resetRands();
funs=bearingCostFunctions('angleSq');
L=10;
NX=10;
TFinal=5;
t=linspace(0,TFinal);
k=1;

X=L*rand(2,NX);
XGoal=L*rand(2,1);
XEval0=[5;5];
dXEval=cnormalize(randn(2,1));
XEval=@(t) XEval0+t*dXEval;

dmax=max(max(evalfun(@(t) computeDi(XEval(t), X, XGoal),t)));
display(dmax)

figure(1)
bearingCostDisplayOld(X,XGoal,funs,L)
xTrack=evalfun(XEval, t);
hold on
plot(xTrack(1,:),xTrack(2,:))
hold off

figure(2)
%check_der(@(t) angle(XEval(t),X,XGoal,dXEval),'function',t)
%plotfun(@(t) computeDi(XEval(t), X, XGoal),t)

plotfun(@(t) dAngle(XEval(t),X,XGoal,dXEval),t,'k')
hold on
plotfun(@(t) dAngleBound(XEval(t),X,XGoal,dXEval,dmax,'lb'),t,'b')
plotfun(@(t) dAngleBound(XEval(t),X,XGoal,dXEval,dmax,'ub'),t,'r')
hold off

function [theta,dtheta]=angle(XEval,X,XGoal,dXEval)
global funs;
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
[theta,gradTheta]=bearingCostGeneral_desiredAngle(YEval,YGoal,funs,nYEval);
dtheta=gradTheta'*dXEval;

function m=computeDi(XEval,X,XGoal)
[~,~,nYEval]=bearingComputeOld(XEval,X,XGoal);
m=(1./nYEval);

function dtheta=dAngle(XEval,X,XGoal,dXEval)
global funs;
[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
[~,gradTheta]=bearingCostGeneral_desiredAngle(YEval,YGoal,funs,nYEval);
dtheta=gradTheta'*dXEval;

function bdtheta=dAngleBound(XEval,X,XGoal,dXEval,dmax,type)
global funs;

[YEval,YGoal,nYEval]=bearingComputeOld(XEval,X,XGoal);
bdtheta=bearingCostGeneral_desiredAngle_derBound(YEval,YGoal,funs,dXEval,dmax,type);

% 
% NY=size(YEval,2);
% c=bearingComputeCosine(YEval,YGoal);
% 
% g=zeros(2,1);
% H=zeros(2,2,NY);
% a=zeros(1,NY);
% 
% for iX=1:NY
%     Yei=YEval(:,iX);
%     Ygi=YGoal(:,iX);
%     ci=c(iX);
%     fci=f(ci);
%     dfci=df(ci);
%     ddfci=ddf(ci);
%     [gradi,H(:,:,iX)]=bearingCostGeneral_terms(Yei,Ygi,fci,dfci,ddfci,ci);
%     g=g+gradi;
% end
% 
% v1=[0 -1; 1 0]*g/(g'*g);
% for iX=1:NY
%     a(iX)=v1'*H(:,:,iX)*dXEval;
% end
% 
% %bdtheta=sum(a./nYEval);
% 
% switch type
%     case 'lb'
%         ktheta=dmax*sum(a(a<0));
%     case 'ub'
%         ktheta=dmax*sum(a(a>0));
% end
% bdtheta=ktheta;
% 
% 
% % sa=sum(a);
% % aBar=a/sa;
% % switch 'lb'
% %     case 'lb'
% %         ktheta=dmax*sum(aBar(aBar<0));
% %     case 'ub'
% %         ktheta=dmax*sum(aBar(aBar>0));
% % end
% % bdtheta=ktheta*sa;
