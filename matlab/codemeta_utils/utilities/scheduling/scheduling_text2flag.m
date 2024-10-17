function flagAvailability=scheduling_text2flag(availability)
flagDebug=false;

[dayList,hourList]=scheduling_dayHours();

reStrDay=reTools_orGroup(dayList);
reStrHour='\d{3,4}';
reStrPeriod=[reStrHour '-' reStrHour];

flagAvailability=false(length(dayList),length(hourList));
for iAvailability=1:length(availability)
    strAvailabilityRow=lower(availability{iAvailability});
    if ~isempty(strAvailabilityRow)
        day=regexpExtract(strAvailabilityRow,reStrDay);
        dayIdx=find(strcmp(day,dayList),1);
        periodList=regexpExtract(strAvailabilityRow,reStrPeriod);
        nbPeriods=length(periodList);
        for iPeriod=1:nbPeriods
            startEndHourStr=pad(regexpExtract(periodList{iPeriod},reStrHour),4,'left','0');
            startHour=str2double(startEndHourStr{1});
            endHour=str2double(startEndHourStr{2});
            startHourIdx=find(hourList>=startHour-eps,1);
            endHourIdx=find(hourList<=endHour+eps,1,'last');
            flagAvailability(dayIdx,startHourIdx:(endHourIdx-1))=true;
        end
        if flagDebug
            disp(strAvailabilityRow)
            disp(hourList(flagAvailability(dayIdx,:))')
        end
    end
end
