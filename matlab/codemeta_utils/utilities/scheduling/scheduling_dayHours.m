function [dayList,hourList,slotMinutes]=scheduling_dayHours()
%The hours in the hour list refer to the *start* of a slot.
dayList={'mon','tue','wed','thu','fri'};%,'sat','sun'};
slotMinutes=30;
hourList=vec((9:17)*100+(0:floor(60/slotMinutes-1))'*slotMinutes)';
