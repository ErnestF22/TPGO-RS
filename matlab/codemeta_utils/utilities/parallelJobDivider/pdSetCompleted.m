function pdSetCompleted(name,idx)
filename=[name '_pdIndicators'];
pdWaitAndLock(filename)
load(filename, 'served', 'completed');

completed(idx)=1;
save(filename, 'served', 'completed');
pdRemoveLock(filename)
