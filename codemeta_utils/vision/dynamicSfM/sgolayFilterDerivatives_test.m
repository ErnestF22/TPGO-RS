function sgolayFilterDerivatives_test

N = 4;                 % Order of polynomial fit
F = 41;                % Window length

dx = 0.2;
xLim = 200;
t = 0:dx:xLim-1;
w=0.04*pi;
A=1;

y = A*sin(w*t);
dy= A*w*cos(w*t);
ddy=-A*w^2*sin(w*t);
%funCheckDerInterpInterp(t,y,dy,t)
%funCheckDerInterpInterp(t,dy,ddy,t)

yNoise=y+0.01*randn(size(y));

[b,g] = sgolay(N,F);   % Calculate S-G coefficients
yFilter=conv(yNoise,g(:,1),'same');
dyFilter=-conv(yNoise,g(:,2),'same')/dx;
ddyFilter=2*conv(yNoise,g(:,3),'same')/(dx^2);

subplot(3,1,1)
plot(t,y,t,yFilter,t,yNoise)
legend('y','yFilter','yNoise')

subplot(3,1,2)
plot(t,dy,t,dyFilter)
legend('dy','dyFilter')

subplot(3,1,3)
plot(t,ddy,t,ddyFilter)
legend('ddy','ddyFilter')
