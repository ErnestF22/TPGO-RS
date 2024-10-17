function rot_interpolationLinear_test
R1=eye(3);
R2=rot([0;0;1]*2*pi/3);
R3=R2*rot([0;1;0]*pi/4);
RReference=@(t) rot_interpolationLinear([0 5 10 15],cat(3,R1,R2,R3,R1),t);
RTraj=RReference(linspace(0,15,100));
plotRotationTrajectory(RTraj)
disp(RTraj)
disp([R1 R2])