function poseEstimationCovariance_comparisonTest
%resetRands()
NTrials=1000;
GEst1=zeros(4,4,NTrials);
Gamma1=zeros(6,6,NTrials);
GEst2=zeros(4,4,NTrials);
Gamma2=zeros(6,6,NTrials);
w=getTextWaitBar(NTrials);
w(0)
NPoints=10;
XCentered=randn(3,NPoints);
XOffset=[0;0;10]*ones(1,NPoints);

GTruth=eye(4);
X1=XCentered+XOffset;
X2=diag([1,1,0.01])*XCentered+XOffset;
x2=projectFromG(GTruth,X2);
x1=projectFromG(GTruth,X1);

for iTrial=1:NTrials
    noise=0.005*randn(size(x1));
    
    x1Noise=x1+noise;
    GEst1(:,:,iTrial)=poseEstimationRefineFromG(GTruth,X1,x1Noise);
    Gamma1(:,:,iTrial)=poseEstimationCovarianceFromG(GEst1(:,:,iTrial),X1,x1Noise);

    x2Noise=x2+noise;
    GEst2(:,:,iTrial)=poseEstimationRefineFromG(GTruth,X2,x2Noise);
    Gamma2(:,:,iTrial)=poseEstimationCovarianceFromG(GEst2(:,:,iTrial),X2,x2Noise);

    w(iTrial)
end

% disp('Estimated poses')
% disp([GEst1 GEst2])

disp('Estimated covariances')
disp(mean(Gamma1,3))
disp(' ')
disp(mean(Gamma2,3))

save([mfilename '_data'])
