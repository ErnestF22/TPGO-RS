function POCposeEstimationPlusRefine
NX=10;
R=rot([pi/8;0;0]);
T=[1;0;0];

X=[randn(2,NX); 5+randn(1,NX)];
x=projectFromRT(R,T,X);

%plot(x(1,:),x(2,:),'*')

[REst,TEst]=poseEstimation(X,[x;ones(1,NX)]);

disp([RT2G(R,T);RT2G(REst,TEst)])
