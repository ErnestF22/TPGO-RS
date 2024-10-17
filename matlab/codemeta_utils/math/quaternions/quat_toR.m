function R=quat_toR(q)
Nq=size(q,2);
if Nq>1
    R=zeros(3,3,Nq);
    for iq=1:Nq
        R(:,:,iq)=quat_toR(q(:,iq));
    end
else
    method='algebraic';
    flagNormalizeQ=true;
    if flagNormalizeQ
        q=cnormalize(q);
    end
    switch method
        case 'angleaxis'
            r=cnormalize(q(2:4,:));
            sthetaHalf=sqrt(sum(q(2:4,:).^2));
            cthetaHalf=q(1,:);
            thetaHalf=atan2(sthetaHalf,cthetaHalf);
            v=r.*([1;1;1]*2*thetaHalf);
            R=rot(v);
        case 'algebraic'
             R=[
                q(1)^2+q(2)^2-q(3)^2-q(4)^2,    -2*q(1)*q(4)+2*q(2)*q(3),       2*q(1)*q(3)+2*q(2)*q(4);
                2*q(1)*q(4) + 2*q(2)*q(3),      q(1)^2-q(2)^2+q(3)^2-q(4)^2,    -2*q(1)*q(2)+2*q(3)*q(4)
                -2*q(1)*q(3)+2*q(2)*q(4),       2*q(1)*q(2) + 2*q(3)*q(4),      q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
    end
end