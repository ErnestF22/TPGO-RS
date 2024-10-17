function [Res] = Colinear(x1,x2,xnode)
tol=1e-6;
temp=0;

a=abs(x2(1)-xnode(1))<tol;
b=abs(xnode(1)-x1(1))<tol;

if a || b
    c=(x1(2)-x2(2))*(xnode(2)-x2(2))>0;
    d=abs(x1(2)-x2(2))>abs(xnode(2)-x2(2));
    if a && b && c && d
        temp=1;
    end
else
    slope1=(x2(2)-xnode(2))/(x2(1)-xnode(1));
    slope2=(xnode(2)-x1(2))/(xnode(1)-x1(1));
    temp = slope2-slope1<tol;
end
Res = temp;
end