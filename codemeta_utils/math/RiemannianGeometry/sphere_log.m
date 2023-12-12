function hi=sphere_log(yi,yj)
%flag to control tricks to improve numerical precision
flagNormalize=false;         %re-normalize result

flagPermute=false;


if(size(yj,2)==1)
    yj=permute(yj,[1 3 2]);
    flagPermute=true;
end
if flagNormalize
    yi=cnormalize(yi);
    yj=cnormalize(yj);
end

[D,N]=size(yj);

hi=zeros(size(yj));
for iN=1:N
    [Q,R]=qr([yi yj(:,iN)],0);
    S=diag(sign(diag(R)));
    Q=Q*S;
    R=S*R;
    hnorm=Q(:,2);
    theta=atan2(R(2,2),R(1,2));
%     [hnorm,r]=qr(sphere_tangentProj(yi,yj(:,iN)),0);
%     if r<0
%         hnorm=-hnorm;
%     end
%     theta=vctAngle(yi,yj(:,iN));
    hi(:,iN)=hnorm*theta;
end

if(flagPermute)
    hi=permute(hi,[1 3 2]);
end
