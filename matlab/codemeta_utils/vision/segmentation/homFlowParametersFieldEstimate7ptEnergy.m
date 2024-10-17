function E=homFlowParametersFieldEstimate7ptEnergy(alphaHat,alpha,x,dx,idxList,kList,varargin)
totalIt=size(alphaHat,3);
lambdaSchedule=zeros(totalIt,3);
lambdaSchedule(:,1)=5;         %cost parameter for data term
lambdaSchedule(:,2)=1e-5;      %cost parameter for soft constraint between alpha and alphaHat
lambdaSchedule(:,3)=1;         %cost parameter for L1 norm
flagDisplay=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'lambdaschedule'
            ivarargin=ivarargin+1;
            lambdaSchedule=varargin{ivarargin};
        case 'display'
            flagDisplay=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

E=zeros(1,totalIt+1);
p=1;
for it=1:totalIt
    lambda=lambdaSchedule(it,:);
    %estimate alphaHat from alpha
    [E(p),E1,E2,E3]=evaluateCost(alphaHat(:,:,it),alpha(:,:,it),x,dx,idxList,kList,lambda);
    fprintfFlag(flagDisplay,'(%d,%d) l = [%.2e %.2e %.2e]; E = %.4f + %.4f + %.4f = %.4f\n',it,p,lambda,E1,E2,E3,E(p));
    p=p+1;
    %estimate alpha from alphaHat at each point
    [E(p),E1,E2,E3]=evaluateCost(alphaHat(:,:,it),alpha(:,:,it+1),x,dx,idxList,kList,lambda);
    fprintfFlag(flagDisplay,'(%d,%d) l = [%.2e %.2e %.2e]; E = %.4f + %.4f + %.4f = %.4f\n',it,p,lambda,E1,E2,E3,E(p));
end

function [E,E1,E2,E3]=evaluateCost(alphaHat,alpha,x,dx,idxList,kList,lambda)
NX=size(x,2);
%add data terms
E1=0;
for iX=1:NX
    ALocal=homFlowParametersEstimate7ptSystemMat(x(:,idxList(iX,1:kList(iX))));
    dxLocal=reshape(dx(:,idxList(iX,1:kList(iX))),[],1);
    E1=E1+lambda(1)*norm(ALocal*alpha(:,iX)-dxLocal)^2;
end
%add constraint term
E2=lambda(2)*norm(alpha-alphaHat,'fro')^2;
%add median regularization term
E3=0;
for iX=1:NX
    E3=E3+lambda(3)*sum(sum(abs(alphaHat(:,idxList(iX,1:kList(iX)))-alphaHat(:,iX)*ones(1,kList(iX)))));
end
E=E1+E2+E3;
