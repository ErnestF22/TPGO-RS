function vH=householderRotation_DDiff(x,k,dx,ddx)
dim=size(x,1);
N=size(x,2);
if N>1
    vH=zeros(dim,N);
    for iN=1:N
        vH(:,iN)=householderRotation_DDiff(x(:,iN),k,dx(:,iN),ddx(:,iN));
    end
else
    I=eye(dim);
    ek=I(:,k);
    xp=cnormalize(x);
    dxp=cnormalizeDiff(x,dx);
    ddxp=cnormalizeDDiff(x,dx,ddx);

    v=xp+ek;
    dv=dxp;
    ddv=ddxp;

    vp=cnormalize(v);
    %dvp=cnormalizeDiff(v,dv);
    ddvp=cnormalizeDDiff(v,dv,ddv);

    vH=-2*hat3(vp)*ddvp;
end
