function [A,b,ref] = fitParallel(n,p,o)
coeff = polyfit([n(1) p(1)],[n(2) p(2)],1);
m = coeff(1);
b0 = coeff(2);
b1 = o(2)-coeff(1)*o(1);
proj_x = ((1-m^2)*o(1)+2*m*o(2)-2*m*b0)/(1+coeff(1)^2);
proj_y = ((m^2-1)*o(2)+2*m*o(1)+2*b0)/(1+coeff(1)^2);

b2 = proj_y-coeff(1)*proj_x;

A = [coeff(1) 1;coeff(1) 1];
b = [b1;b2];
ref = [proj_x,proj_y];
end
    