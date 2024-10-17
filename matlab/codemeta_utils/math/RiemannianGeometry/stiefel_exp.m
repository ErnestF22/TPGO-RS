%Compute exponential map for a vector H in the tangent space of the point Y
%on the Stiefel manifold (see Edelman et. al. for referring equation numbers)
function [Yt,Ht,A,Q,R,MN,MNbig]=stiefel_exp(Y,H)

N2=size(H,3);
if(N2>1)
    %recursive call if H contains multiple tangent vectors
    Yt=zeros([size(Y) N2]);
    if(nargout>1)
        Ht=zeros([size(Y) N2]);
    end
    for(iN2=1:N2)
        if(nargout==1)
            Yt(:,:,iN2)=stiefel_exp(Y,H(:,:,iN2));
        else
            [Yt(:,:,iN2),Ht(:,:,iN2)]=stiefel_exp(Y,H(:,:,iN2));
        end
    end
else
%     %make sure H is in the tangent space
%     H=stiefel_tangentProj(Y,H);
    
    %compute the exponential
    [~,p]=size(Y);

    A=[Y'*H];
    A=(A-A')/2; %make sure it is skew-symmetric
    [Q,R]=qr(H - Y*A,0);            %H=Y*A+Q*R
    MN=rot_exp(eye(2*p), [A, -R'; R, zeros(p)]);
    MNbig=MN;
    MN=MN(:,1:p);
    Yt=Y*MN(1:p,:)+Q*MN(p+1:2*p,:); %Geodesic from (2.45)
    if(nargout>1)
        Ht=H*MN(1:p,:)-Y*(R'*MN(p+1:2*p,:)); %Geodesic direction from (3.3) at Yt (from parallel transport)
    end
end

