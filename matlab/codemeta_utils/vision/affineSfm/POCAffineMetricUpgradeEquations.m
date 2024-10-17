function POCAffineMetricUpgradeEquations
M=sym('M',[2 3]);
%g=[G11 G12 G13 G22 G23 G33]
g=sym('g',[6 1]);
G=[g(1) g(2) g(3); g(2) g(4) g(5); g(3) g(5) g(6)];
disp(G-G.')

E=M*G*M.';

[A1,T1]=coeffs(E(1,1),g);
b1=1;
assert(all((T1-g.')==0))

[A2,T2]=coeffs(E(1,2),g);
b2=0;
assert(all((T2-g.')==0))

[A3,T3]=coeffs(E(2,2),g);
b3=1;
assert(all((T3-g.')==0))

A=[A1;A2;A3];
b=[1;0;1];
fileNameFunction='affineMetricUpgradeCoefficients';
matlabFunction(A,'File','affineMetricUpgradeCoefficients','Vars',{M},'Optimize',false);
disp(['Code for computing matrix of coefficients saved to ' fileNameFunction])

