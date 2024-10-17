function scheduling_input()

name=input('Name: ','s');
reqTime=input('Required time (in minutes, multiple of 15): ','s');
inStr='bu';
availability={};
disp('Availability:')
while ~isempty(inStr)
    inStr=input('','s');
    availability{end+1}=inStr;
end

save(fullfile('schedules',name),'name','reqTime','availability')


