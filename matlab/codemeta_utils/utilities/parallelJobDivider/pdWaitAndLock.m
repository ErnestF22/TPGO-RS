function pdWaitAndLock(name)
lockName=[name '_lock.mat'];
t=clock;
tic
nextWarningIncrement=10;
nextWarningTime=nextWarningIncrement;

%poll for the lock
while exist(lockName,'file')
    tWaiting=toc;
    if tWaiting>nextWarningTime
        disp(['Waiting for lock ' lockName ' for more than ' num2str(tWaiting) ' seconds'])
        nextWarningTime=nextWarningTime+nextWarningIncrement;
    end        
    pause(0.01)
end

save(lockName,'t')
