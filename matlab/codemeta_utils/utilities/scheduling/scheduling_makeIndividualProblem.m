function problem=scheduling_makeIndividualProblem(flag,slotLength)
[dayList,hourList]=scheduling_dayHours();
nbDays=length(dayList);
nbHours=length(hourList);

[~, basisWeekDiff]=scheduling_basisWeek();

nbTotalSlots=nbDays*nbHours;
u=ones(nbTotalSlots,1);
z=zeros(nbTotalSlots,1);

%disp(norm(flag-full(basis2flag(flag2basis(flag),nbDays)),'fro'))

AeqLength=ones(1,nbTotalSlots)*basisWeekDiff;
beqLength=slotLength;
AeqUniqueStart=sparse([u;z])';
beqUniqueStart=1;
AeqUniqueEnd=sparse([z;u])';
beqUniqueEnd=1;
Aeq=[AeqLength;AeqUniqueStart;AeqUniqueEnd];
beq=[beqLength;beqUniqueStart;beqUniqueEnd];

AineqAvailability=basisWeekDiff;
bineqAvailability=scheduling_flag2basis(flag);
AineqPositiveOrder=-basisWeekDiff;
bineqPositiveOrder=zeros(nbTotalSlots,1);

Aineq=[AineqAvailability;AineqPositiveOrder];
bineq=[bineqAvailability;bineqPositiveOrder];

problem.solver='intlinprog';
problem.options=optimoptions(problem.solver);
problem.f=basisWeekDiff'*ones(nbTotalSlots,1);
problem.intcon=1:(2*nbTotalSlots);
problem.lb=zeros(2*nbTotalSlots,1);
problem.ub=ones(2*nbTotalSlots,1);
problem.Aeq=Aeq;
problem.beq=beq;
problem.Aineq=Aineq;
problem.bineq=bineq;
