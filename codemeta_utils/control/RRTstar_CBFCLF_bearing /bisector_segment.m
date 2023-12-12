% input: two points a in R2 and b in R2
% Computes the bisector of a segment line between point a and b
% output: A*[x;-y]+B=0
function [A,B] = bisector_segment(a,b)
c = (a+b)/2;
coeff = polyfit([a(1) b(1)],[a(2) b(2)],1);
if abs(coeff(1))>=10^5 %vertical line
    A=[0 1];
    B=c(2);
else
    if abs(coeff(1))<10^-5%horizontal line
        A=[1 0];
        B=c(1);
    else
        A=[-1/coeff(1) 1];
        B=c(2)-(-1/coeff(1))*c(1);
    end   
end
end