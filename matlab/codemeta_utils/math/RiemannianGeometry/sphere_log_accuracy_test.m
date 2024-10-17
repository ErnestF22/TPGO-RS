function sphere_log_accuracy_test
resetRands();
n=7;
NTrials=1000;
y1=cnormalize(randn(n,1));

%sigmanoise=1e-12;
% err=NaN(NTrials,1);
% for iTrial=1:NTrials
%     H=sphere_tangentProj(y1,randn(n,1));
%     H=H/sqrt(sphere_metric(y1,H,H));
%     y2=cnormalize(sphere_exp(y1,sigmanoise*randn*H));
% %     y2=cnormalize(randn(n,1));
%     err(iTrial)=sum((y2-sphere_exp(y1,sphere_log(y1,y2))).^2);
% end

y2=cnormalize(randn(n,NTrials));
err=sqrt(sum((y2-sphere_exp(y1,sphere_log(y1,y2))).^2));

plot(err)
disp(['Mean: ' num2str(mean(err))])
disp(['Max: ' num2str(max(err))])
disp(['Std:  ' num2str(std(err))])