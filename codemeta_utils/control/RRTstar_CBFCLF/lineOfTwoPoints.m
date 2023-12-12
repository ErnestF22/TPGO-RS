function [A,B]=lineOfTwoPoints(X1,X2)
%finding the prependecular line with passes through X1 and X2:
coeff = polyfit([X1(1) X2(1)], [X1(2) X2(2)], 1);
xMid=(X1(1)+X2(1))/2;
yMid=(X1(2)+X2(2))/2;
if abs(coeff(1))<=10^-5 %vertical line
    A=[0 1];
    B=yMid;
else
    if abs(coeff(1))>10^5%horizontal line
        A=[1 0];
        B=xMid;
    else
        A=[coeff(1) 1];
        B=coeff(2);
    end   
end
end
