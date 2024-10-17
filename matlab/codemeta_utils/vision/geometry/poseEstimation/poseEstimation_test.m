function poseEstimation_test

resetRands();

Ntrials=5;
Nposes=3;

Lnoises=[0 0.1 0.5 1 2 3];
sigmanoise=0.01;

Txmax=10;
Tymax=10;
Tzmax=100;
Tzmin=5;

%3D Model
% X=[ -1  0  0;
%      1  0  0;
%      0  1  0;
%      0 -1  0;
%      0  0 -1;
%      0  0  1]';

X=randn(3,15);


zmaxmin=max(X(3,:))-min(X(3,:));

%n=size(X,2);    %number of points in the model

%randn('state',0)
%rand('state',0)

randn('state',0)
rand('state',0)

for(ipose=1:Nposes)
    fprintf('Pose %d\n',ipose)
    %Pose ground truth
    R_truth=rot(randn(3,1));
    T_truth=[2*Txmax*(rand-0.5);2*Tymax*(rand-0.5);Tzmin*zmaxmin+rand*(Tzmax-Tzmin)*zmaxmin];
    
    allPoses{ipose}.R=R_truth;
    allPoses{ipose}.T=T_truth;

    %Transform and project
    p=projectFromRT(R_truth,T_truth,X,'poses');

    sigmap=sqrt(var(p',1));
    sigmap=mean(sigmap(1:2));

    for itrial=1:Ntrials
        fprintf('Trial %d\n',itrial)
        noise=randn(2,size(p,2));
        for inoise=1:length(Lnoises)
            sigmaeff=Lnoises(inoise)*sigmanoise*sigmap;

            p_noise=p;
            p_noise(1:2,:)=p(1:2,:)+max(min(sigmaeff*noise,3*sigmaeff),-3*sigmaeff);

            for alg=1:3
                switch(alg)
                    case 1
                        GEst=poseEstimationQuadratic(X,p_noise,'poses');
                    case 2
                        GEst=poseEstimationDLT(X,p_noise);
                    case 3
                        GEst=poseEstimation(X,p_noise,'poses');
                    case 4
                        GEst=poseEstimationDLT(X,p_noise,'secondorder');
                end
                [R_est,T_est]=G2RT(GEst);
                allT_est{alg,inoise,itrial,ipose}=T_est;
                allR_est{alg,inoise,itrial,ipose}=R_est;
                roterror(alg,inoise,itrial,ipose)=rot_dist(R_truth,R_est);
                translerror(alg,inoise,itrial,ipose)=norm(T_truth-T_est)/(norm(T_truth)+norm(T_est))*2;
            end
        end
    end
end
figure
plot(Lnoises*sigmanoise,mean(mean(roterror,4),3)')
title('Rotation error');
legend('Ansar','Hartley','Ansar+Refine')

figure
plot(Lnoises*sigmanoise,mean(mean(translerror,4),3)')
title('Relative translation error');
legend('Ansar','Hartley','Ansar+Refine')

save([mfilename '_Nposes' num2str(Nposes) '_Ntrials' num2str(Ntrials)])
