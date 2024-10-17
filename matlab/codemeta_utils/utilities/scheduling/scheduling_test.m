function scheduling_test
[data,dataOrder]=scheduling_loadAvailability();
data=scheduling_addIndividualProblems(data);
data.Tron.f=scheduling_baseCost(data);

problem=scheduling_makeGlobalProblem(data,dataOrder);
xSol=intlinprog(problem);

save([mfilename '_testData'])

nbDays=length(scheduling_dayHours());
[~, basisWeekDiff]=scheduling_basisWeek();
idxStudent=reshape(1:problem.szDecision,problem.szDecisionStudent,problem.nbStudents);
slotStr=cell(1,problem.nbStudents-1);
nameList=dataOrder(2:end);
for iStudent=1:problem.nbStudents
    name=nameList{iStudent};
    xStudent=xSol(idxStudent(:,iStudent));
    basis=basisWeekDiff*xStudent;
    flag=scheduling_basis2flag(basis,nbDays);
    slotStr{iStudent}=char(scheduling_flag2text(flag));
%     disp(['# ' name])
%     disp(slotStr{iStudent})
%     disp(strjoin(data.(name).availability,','))
end

[slotStrSorted,idxSlotStrSorted]=sort(slotStr);
nameListSorted=nameList(idxSlotStrSorted);
disp([char(nameListSorted) ' '*ones(problem.nbStudents,1) char(slotStrSorted)])



function [problem,szDecisionStudent]=scheduling_makeGlobalProblem(data,dataOrder)
nbStudents=length(dataOrder)-1;
[~, basisWeekDiff]=scheduling_basisWeek();
szDecisionStudent=size(basisWeekDiff,2);
szDecision=szDecisionStudent*nbStudents;
nameBase=dataOrder{1};
AeqList=cell(1,nbStudents);
beqList=cell(1,nbStudents);
AineqList=cell(1,nbStudents);
bineqList=cell(1,nbStudents);

for iStudent=1:nbStudents
    name=dataOrder{iStudent+1};
    AeqList{iStudent}=data.(name).problem.Aeq;
    beqList{iStudent}=data.(name).problem.beq;
    AineqList{iStudent}=data.(name).problem.Aineq;
    bineqList{iStudent}=data.(name).problem.bineq;
end

AineqBase=repmat(basisWeekDiff,1,nbStudents);
bineqBase=scheduling_flag2basis(data.(nameBase).flag);
if isfield(data.(nameBase),'f')
    f=AineqBase'*scheduling_flag2basis(data.(nameBase).f);
else
    f=ones(szDecision,1);
end

problem.Aeq=blkdiag(AeqList{:});
problem.beq=cat(1,beqList{:});
problem.Aineq=[AineqBase; blkdiag(AineqList{:})];
problem.bineq=cat(1,bineqBase,bineqList{:});
problem.f=f;
problem.intcon=1:(szDecision);
problem.lb=zeros(szDecision,1);
problem.ub=ones(szDecision,1);
problem.solver='intlinprog';
problem.options=optimoptions(problem.solver);

problem.nbStudents=nbStudents;
problem.szDecisionStudent=szDecisionStudent;
problem.szDecision=szDecision;

function data=scheduling_addIndividualProblems(data)
nameList=fieldnames(data);
for iName=1:length(nameList)
    name=nameList{iName};
    [data.(name).problem]=scheduling_makeIndividualProblem(data.(name).flag,data.(name).reqTime);
end




