function scheduling_makeIndividualProblem_test
data=scheduling_loadAvailability();
disp(data)

[data.Cheng.problem]=scheduling_makeIndividualProblem(data.Cheng.flag,3);

xSol=intlinprog(data.Cheng.problem);

nbDays=length(scheduling_dayHours());
[~, basisWeekDiff]=scheduling_basisWeek();
disp(full(scheduling_basis2flag(basisWeekDiff*xSol,nbDays)>0.1))
