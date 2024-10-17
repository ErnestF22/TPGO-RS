function availability=scheduling_flag2text(flag)
[dayList,hourList]=scheduling_dayHours();

flagDiff=diff(flag,[],2);
availability={};
for iDay=1:length(dayList)
    startSlot=find(flagDiff(iDay,:)==1)+1;
    endSlot=find(flagDiff(iDay,:)==-1)+1;
    
    %print day only if there is at least one slot
    if ~isempty(startSlot) || ~isempty(endSlot)
        %check if the start/end fall at any of the edges of the day
        if length(startSlot)<length(endSlot)
            startSlot=[1 startSlot];
        end
        if length(endSlot)<length(startSlot)
            endSlot=[endSlot length(hourList)];
        end
        
        availabilityDay=sprintf('%s',dayList{iDay});
        for iPeriod=1:length(startSlot)
            availabilityDay=[availabilityDay sprintf(', %04d-%04d', hourList(startSlot(iPeriod)), hourList(endSlot(iPeriod)))];
        end
        availability{end+1}=availabilityDay;
    end
end
