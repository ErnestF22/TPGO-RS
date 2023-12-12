function x = find_intersecition(p1,q1,p2,q2)
coeff_1 = polyfit([p1(1) q1(1)],[p1(2) q1(2)],1);
coeff_2 = polyfit([p2(1) q2(1)],[p2(2) q2(2)],1);
x_1 = (coeff_2(2)-coeff_1(2))/(coeff_1(1)-coeff_2(1));
x_2 = coeff_1(1)*x_1+coeff_1(2);
x = [x_1; x_2];
end