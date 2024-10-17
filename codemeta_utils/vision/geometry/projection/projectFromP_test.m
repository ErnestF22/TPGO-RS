function projectFromP_test
s=load('triangulate_test_dataset_data.mat');
P=RTK2P(s.R(:,:,2),s.T(:,2),s.K);

dX=randn(3,1);
w=randn(3,1);
X=@(t) [0;0;10]+t*dX;

figure(1)
check_der(@(t) evaluateX(P,X,dX,t),'function');
 
figure(2)
check_der(@(t) evaluatedX(P,X,dX,w,t),'function');

function [x,dx]=evaluateX(P,X,dX,t)
[x,Jx]=projectFromP(P,X(t));
dx=Jx*dX;

function [dx,ddx]=evaluatedX(P,X,dX,w,t)
[~,Jx,Hx]=projectFromP(P,X(t));
dx=Jx*w;
ddx=zeros(2,1);
for k=1:2
    ddx(k)=dX'*Hx(:,:,k)*w;
end
