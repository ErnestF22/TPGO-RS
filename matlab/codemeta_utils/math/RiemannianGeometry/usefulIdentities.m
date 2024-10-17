function usefulIdentities
%Numerical verification of a few identities

%hatProperties()
figure(1)
compareFunsTheta('(cos(theta)+1)/sin(theta)','cot(theta/2)');
figure(2)
compareFunsTheta('(cos(theta)-1)/sin(theta)','-tan(theta/2)');
figure(3)
%compareFunsTheta('sqrt((theta/2)^2+(theta/2)^2*cot(theta/2)^2)','csc(theta)');
%compareFunsTheta('sqrt((theta/2)^2+(theta/2)^2*cot(theta/2)^2)','abs(theta/2)*sqrt(1+cot(theta/2)^2)')
compareFunsTheta('sqrt(1+cot(theta/2)^2)','csc(theta/2)')
figure(4)
compareFunsTheta('sin(2*theta)/sin(theta)','2*cos(theta)');

function compareFunsTheta(sf1,sf2)
t=linspace(0, 0.99*pi,200);
f1=@(theta) eval(sf1);
f2=@(theta) eval(sf2);
plotfun(f1,t,'gx')
hold on
plotfun(f2,t,'b')
hold off
legend(sf1,sf2)

N=4; M=2;
IN=eye(N); IM=eye(M); P=randn(M,N); Q=randn(N,M); 
disp('P*inv(I+Q*P)-inv(I+P*Q)*P')
disp(P*inv(IN+Q*P)-inv(IM+P*Q)*P)
