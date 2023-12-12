function SVD_test(angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x=[cos(angle); sin(angle)];
dx=[-0.1*sin(angle); 0.1*cos(angle)];
landmarks=[-6.3267 -1.4564 6.4309 -9.4093 -5.0582 -1.3374;...
            -5.6740 9.4115 -2.6130 -6.1627 1.3438 2.2212];
[y,ny]=bearingCompute(x,landmarks);
dy=bearingComputeDerivative(dx,y,ny);
landmarks=bearingComputeSVDEstimate(y,dy,dx)
m1=y(2,1)/y(1,1);
m2=y(2,2)/y(1,2);
X(1)=(m1*landmarks(1,1)-m2*landmarks(1,2)-landmarks(2,1)+landmarks(2,2))/(m1-m2);
X(2)=m1*(X(1)-landmarks(1,1))+landmarks(2,1);



