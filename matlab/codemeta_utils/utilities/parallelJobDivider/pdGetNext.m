function idx=pdGetNext(name)
filename=[name '_pdIndicators'];

pdWaitAndLock(filename)
load(filename, 'served', 'completed');

idx=find(~served,1,'first');

if ~isempty(idx)
    served(idx)=1;
    save(filename, 'served', 'completed');
end
pdRemoveLock(filename)
