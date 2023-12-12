function scheduling_loadAvailability_test()
data=scheduling_loadAvailability();

nameList=fieldnames(data);
for iName=1:length(nameList)
    name=nameList{iName};
    disp(['# ' name])
    disp(scheduling_flag2text(data.(name).flag)')
    s=load(['schedules/' name '.mat']);
    disp(s.availability')
end