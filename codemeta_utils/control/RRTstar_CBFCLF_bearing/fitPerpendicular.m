function [A,B] = fitPerpendicular(n,p)
coeff = polyfit([n(1) p(1)],[n(2) p(2)],1);
A=[-1/coeff(1) 1];
B=n(2)-(-1/coeff(1))*n(1);
end