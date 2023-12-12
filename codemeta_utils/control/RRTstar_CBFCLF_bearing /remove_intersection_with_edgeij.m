function [A,b] = remove_intersection_with_edgeij(A,b,n,p)
coeff =  polyfit([n(1) p(1)],[n(2) p(2)],1);
m0 = coeff(1);
b0 = coeff(2);
idx = [];
for i=1:size(A,1)
    m1 = -A(i,1)/A(i,2);
    b1 = b(i)/A(i,2);
    x = (b1-b0)/(m0-m1);
    r = (x-p(1))/(n(1)-p(1));
    if (r>=0) && (r<=1)
        idx = [idx i];
    end
end
A(idx,:)=[];
b(idx,:)=[];
end