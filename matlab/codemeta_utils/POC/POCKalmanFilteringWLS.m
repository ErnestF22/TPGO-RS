function POCKalmanFilteringWLS

x0=rand;
s0=rand;
x1=rand;
s1=rand;

S=s0+s1;
xhatKF=(1-s0/S)*x0+s0/S*x1;
xhatWLS=(x0/s0+x1/s1)/(1/s0+1/s1);

disp([xhatKF xhatWLS])
