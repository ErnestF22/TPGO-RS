function [data,dataOrder]=scheduling_loadAvailability()
data=[];
[~,~,slotMinutes]=scheduling_dayHours();

listing=dir(fullfile('schedules','*.mat'));
nameList=strrep({listing.name},'.mat','');

%load each available schedule
for iName=1:length(nameList)
    name=nameList{iName};
    s=load(fullfile('schedules',[name '.mat']));
    availability=strrep(s.availability,':','');
    data.(name).availability=availability;
    data.(name).flag=scheduling_text2flag(availability);
    data.(name).reqTime=floor(str2double(s.reqTime)/slotMinutes);
end

%order
dataOrder=nameList;
idxBase=find(strcmpi(nameList,'Tron'),1);
dataOrder([1 idxBase])=dataOrder([idxBase 1]);

for iOrder=1:length(dataOrder)
    data.(dataOrder{iOrder}).order=iOrder;
end
