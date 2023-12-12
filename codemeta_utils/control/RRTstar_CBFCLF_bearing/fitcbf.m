function [A,b] = fitcbf(n,p,o)
coeff = polyfit([n(1) p(1)],[n(2) p(2)],1);
m = coeff(1);
b0 = coeff(2);
proj_x = ((1-m^2)*o(1)+2*m*o(2)-2*m*b0)/(1+coeff(1)^2);
proj_y = ((m^2-1)*o(2)+2*m*o(1)+2*b0)/(1+coeff(1)^2);

x = [proj_x;proj_y];

coeff = polyfit([x(1) p(1)],[x(2) p(2)],1);
A1 = [coeff(1) 1];
b1 = x(2)-coeff(1)*x(1);


coeff = polyfit([o(1) p(1)],[o(2) p(2)],1);
A2 = [coeff(1) 1];
b2 = o(2)-coeff(1)*o(1);

A = [A1;A2];
b = [b1;b2];
end
    