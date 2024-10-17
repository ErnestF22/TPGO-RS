function [Yt,U,S,V]=grassman_exp(Y,H)
%recursive call if H contains multiple tangent vectors
N2=size(H,3);
if(N2>1)
    Yt=zeros([size(Y) N2]);
    for iN2=1:N2
        Yt(:,:,iN2)=grassman_exp(Y,H(:,:,iN2));
    end
else
    [U,S,V]=svd(H,'econ');
    Yt=[Y*V U]*[diag(cos(diag(S)));diag(sin(diag(S)))]*V';
    [Yt,void]=qr(Yt,0);
end
