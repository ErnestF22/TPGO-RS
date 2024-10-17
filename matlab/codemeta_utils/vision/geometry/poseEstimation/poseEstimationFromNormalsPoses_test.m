function poseEstimationFromNormalsPoses_test
load poseEstimationFromNormalsPoses_dataset

e3=[0;0;1];
%check constraints
disp(max(abs(Rn(:,:,1)*e3-RTruth*n(:,1))))
disp(Rn(:,3,1)'*TTruth-Rn(:,3,1)'*Tn(:,1)+d(1))

[R,T]=poseEstimationFromNormalsPoses(n,d,Rn,Tn);

disp([R RTruth])
disp([T TTruth])

