function angleDiff_test
theta1=pi;
theta2=linspace(-4*pi,4*pi,600);

thetaDiff=angleDiff(theta1,theta2);

plot(theta2/pi,thetaDiff/pi)
grid on