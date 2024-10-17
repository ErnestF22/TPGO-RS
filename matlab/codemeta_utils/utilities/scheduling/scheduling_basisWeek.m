function [basisWeek, basisWeekDiff]=scheduling_basisWeek()
%basisWeekDiff is the matrix going from the start/end decision variables to the flags (vectorized).
[dayList,hourList]=scheduling_dayHours();
nbDays=length(dayList);
nbHours=length(hourList);

basisDay=sparse(triu(ones(nbHours)));
basisWeek=kron(speye(nbDays),basisDay);
% for iBase=1:48
%     disp(full(basis2flag(basisWeek(:,iBase),nbDays)))
%     pause
% end

%matrix to go from (day,slot) to basis
idxDaySlot=reshape(1:nbDays*nbHours,nbHours,nbDays)';
%disp(full(basis2flag(basisWeek(:,idxDaySlot(2,3)),nbDays)))

day=randi(nbDays);
slotStart=randi(nbHours);
% xStart=zeros(nbDays*nbSlots,1);
% xStart(idxDaySlot(day,slotStart))=1;
% xEnd=zeros(nbDays*nbSlots,1);
% xEnd(idxDaySlot(day,slotStart+1))=1;
% fprintf('day: %d, slotStart: %d\n',day,slotStart)
% disp(full(basis2flag([-basisWeek basisWeek]*[xStart;xEnd],nbDays)))
basisWeekDiff=[-basisWeek basisWeek];
