function [q,d]=rot2quat(R)
for(iR=1:size(R,3))
    Qxx=R(1,1,iR); Qxy=R(1,2,iR); Qxz=R(1,3,iR);
    Qyx=R(2,1,iR); Qyy=R(2,2,iR); Qyz=R(2,3,iR);
    Qzx=R(3,1,iR); Qzy=R(3,2,iR); Qzz=R(3,3,iR);
    if(Qxx>Qyy && Qxx>Qzz)
        r = sqrt(1+Qxx-Qyy-Qzz);
        s = 0.5/r;
        w(1,iR) = (Qzy-Qyz)*s;
        x(1,iR) = 0.5*r;
        y(1,iR) = (Qxy+Qyx)*s;
        z(1,iR) = (Qzx+Qxz)*s;
        d(1,iR) = 1;
    else
        if(Qyy>Qxx && Qyy>Qzz)
            r = sqrt(1+Qyy-Qxx-Qzz);
            s = 0.5/r;
            w(1,iR) = (Qxz-Qzx)*s;
            x(1,iR) = (Qxy+Qyx)*s;
            y(1,iR) = 0.5*r;
            z(1,iR) = (Qzy+Qyz)*s;
            d(1,iR) = 2;
        else
            r = sqrt(1+Qzz-Qxx-Qyy);
            s = 0.5/r;
            w(1,iR) = (Qyx-Qxy)*s;
            x(1,iR) = (Qxz+Qzx)*s;
            y(1,iR) = (Qzy+Qyz)*s;
            z(1,iR) = 0.5*r;
            d(1,iR) = 3;
        end
    end
end

q=[w;x;y;z];
q=q.*repmat(sign(q(1,:)),4,1);

