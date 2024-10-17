function [EE,C,membership]=Factorgraph_removeTriangle2(EE,C,membership)
NPins=max(EE(:,1));
%disp('Triangle found');
CycleComps=EE(C,2);
NewComp=min(CycleComps);

for idxCycleComps=1:length(CycleComps)    
    EE(EE(:)==CycleComps(idxCycleComps))=NewComp;
    membership(membership(:)==CycleComps(idxCycleComps)-NPins)=NewComp-NPins;
end

EE=unique(EE,'rows');
end