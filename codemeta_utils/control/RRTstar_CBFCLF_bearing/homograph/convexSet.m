function [A,b]=convexSet(y)
A=[];
b=[];
for i=1:size(y,2)-1
    p1 = y(:,i);
    p2 = y(:,i+1);
    coeff = polyfit([p1(1), p2(1)], [p1(2), p2(2)], 1);
    A = [A,coeff(1)];
    b = [b,coeff(2)];
end
p1 = y(:,1);
p2 = y(:,end);
coeff = polyfit([p1(1), p2(1)], [p1(2), p2(2)], 1);
A = [A,coeff(1)]';
b = [b,coeff(2)]';
A = [A ones(size(A,1),1)];
% A = [a,1] where y = ax+b 
end