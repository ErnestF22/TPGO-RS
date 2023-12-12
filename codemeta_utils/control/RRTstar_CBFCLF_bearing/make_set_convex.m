function [Ax,bx]=make_set_convex(A,b,p)
Ax=A;
bx=b;
for i=1:size(A,1)
    if A(i,2)==0
        d=sign(p(1)-(b(i)));
        if d>0
            Ax(i,1)=-Ax(i,1);
            bx(i)=-bx(i);
        end
    else
        d=sign(p(2)-(p(1)*A(i,1)+b(i)));
        if d<0
            Ax(i,1)=-Ax(i,1);
        else
            Ax(i,2)=-Ax(i,2);
            bx(i)=-bx(i);
        end
    end
end
end