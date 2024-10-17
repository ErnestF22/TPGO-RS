function projectFromRT_test
resetRands(3)
X=randn(3,1);
methodPoses='reference';
%methodPoses='pose';
[Rt,vt,R0,v0,vVec]=rot_randGeodFun(eye(3));
[Tt,dTt,T0,dT]=real_randGeodFun(randn(3,1));
wFix=rand(6,1);

figure(1)
check_der(@(t) evaluateFandDF(X,Rt,Tt,vVec,dT,methodPoses,t),'function')

figure(2)
check_der(@(t) evaluateDFandDdF(X,Rt,Tt,vVec,dT,methodPoses,wFix,t),'function')

function [f,df]=evaluateFandDF(X,Rt,Tt,vVec,dT,methodPoses,t)
[f,Jf]=projectFromRT(Rt(t),Tt(t),X,'methodAbsolutePoses',methodPoses);
% R=Rt(t);
% T=Tt(t);
% [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,'methodAbsolutePoses',methodPoses);
% xTransformed=XTransformed(1:2,:)./repmat(XTransformed(3,:),2,1);
% e3=[0;0;1];
% Jf=1/XTransformed(3)^2*[eye(2) zeros(2,1)]*(e3'*XTransformed*eye(3)-XTransformed*e3')*JXTransformed;
df=Jf*[vVec;dT];

function [df,ddf]=evaluateDFandDdF(X,Rt,Tt,vVec,dT,methodPoses,wFix,t)
[~,Jf,Hf]=projectFromRT(Rt(t),Tt(t),X,'methodAbsolutePoses',methodPoses);
w=[vVec;dT];
df=Jf*wFix;
ddf=zeros(2,1);
for k=1:2
    ddf(k)=w'*Hf(:,:,k)*wFix;
end

% w=[vVec;dT];
% R=Rt(t);
% T=Tt(t);
% [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,'methodAbsolutePoses',methodPoses);
% xTransformed=XTransformed(1:2,:)./repmat(XTransformed(3,:),2,1);
% Xi=XTransformed;
% Ji=JXTransformed;
% Hi=HXTransformed;
% e3Xi=XTransformed(3);
% P=[eye(2) zeros(2,1)];
% E=eye(3);
% e3=E(:,3);
% f1=1/e3Xi^2;
% f2=P*(e3'*Xi*eye(3)-Xi*e3');
% f3=Ji*wFix;
% df=f1*f2*f3;
% 
% 
% df1f2f3=zeros(2,1);
% ddf2f3=zeros(2,1);
% ddf3a=zeros(2,1);
% ddf=zeros(2,1);
% for k=1:2
%     ek=E(:,k);
%     df1f2f3(k)=w'*((-2/e3Xi^3)*Ji'*e3*f2(k,:)*Ji)*wFix;
%     ddf2f3(k)=w'*f1*Ji'*(e3*ek'-ek*e3')*Ji*wFix;
%     f2H=zeros(6,6);
%     for l=1:3
%         f2H=f2H+f2(k,l)*Hi(:,:,l);
%     end
%     ddf3a(k)=w'*f1*f2H*wFix;
%     ddf(k)=df1f2f3(k)+ddf2f3(k)+ddf3a(k);
% end
% 
% % df2=(e1'*eb*e3'*Ji*w-e3'*eb*e1'*Ji*w);
% % ddf=df2;
% % %df=f3;
% % ddf3=zeros(3,1);
% % for l=1:3
% %     ddf3(l)=w'*Hi(:,:,l)*wFix;
% % end
% 
% 
% %ddf=df1+ddf2f3+f1*f2*ddf3;
% ddf=df1f2f3+ddf2f3+ddf3a;
% 
% % Jf=([eye(2) -xTransformed(:)]*Ji(:,:))./e3Xi;
% % df=Jf*w;
% 
% % m1=1/e3Xi^2;
% % m2=[eye(2) zeros(2,1)]*(e3'*Xi*eye(3)-Xi*e3');
% % m3=Ji*wFix;
% % df=m1*m2*m3;
% % 
% % ddf1=-2/e3Xi^3*e3'*Ji*w*m2*m3;
% % ddf2=zeros(2,1);
% % for k=1:2
% %     ek=E(:,k);
% %     ddf2(k)=m1*(ek'*m3*e3'-e3'*m3*ek')*Ji*w;
% % end
% % ddf3=zeros(2,1);
% % for k=1:2
% %     ddf3(k)=w'*HXTransformed(:,:,k)*wFix;
% % end
% % ddf=ddf1;%+ddf2+ddf3;

