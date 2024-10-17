function essential_distMinAngle_speedTest
fileName=[mfilename '_data.mat'];
if ~exist(fileName,'file')
    NTrials=20000;
    Q1=zeros(6,3,NTrials);
    Q2=zeros(6,3,NTrials);
    for iTrial=1:NTrials
        Q10=rot_randn([],[],2);
        Q20=rot_randn([],[],2);
        Q1(:,:,iTrial)=[Q10(:,:,1);Q10(:,:,2)];
        Q2(:,:,iTrial)=[Q20(:,:,1);Q20(:,:,2)];
    end
    save(fileName)
else
    load(fileName)
end

t0=cputime;
for iTrial=1:NTrials
    [tMin,fMin]=essential_distMinAngle(Q1(:,:,iTrial),Q2(:,:,iTrial));
end
tf=cputime;
disp('Total time (s)')
disp(tf-t0)
disp('Time per instance (ms)')
disp((tf-t0)/NTrials*1000)
