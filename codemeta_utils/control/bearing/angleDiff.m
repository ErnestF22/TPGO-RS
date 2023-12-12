function thetaDiff=angleDiff(theta1,theta2)
thetaDiff=mod(theta1-theta2+pi,2*pi)-pi;
