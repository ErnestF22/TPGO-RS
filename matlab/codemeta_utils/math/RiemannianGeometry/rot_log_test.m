randn('state',0)
n=7;
Ntrials=5000;
e=zeros(1,Ntrials);
tic
for trial=1:Ntrials;
    R1=orth(randn(n));
    %R1=eye(n);
    A=rot_tangentProj(R1,randn(n));
    A=A/sqrt(rot_metric(R1,A,A))*0.1*pi;
    e(trial)=norm(A-rot_log(R1,rot_exp(R1,A)),'fro');
end
toc
disp(['mean ' num2str(mean(e)) ' std ' num2str(std(e))])
figure(1)
hist(e)
%figure(1)
%plot(e)
