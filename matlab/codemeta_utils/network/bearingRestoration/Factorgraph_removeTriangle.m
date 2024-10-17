function [EE,C,membership]=Factorgraph_removeTriangle(EE,C,membership)
NPins=max(EE(:,1));
while sum(C(1,:))==6 || sum(C(1,:))==4
    %disp('Triangle found');
    CycleComps=EE(logical(C(1,:)),2);
    NewComp=min(CycleComps);

    for idxCycleComps=1:length(CycleComps)    
        EE(EE(:)==CycleComps(idxCycleComps))=NewComp;
        membership(membership(:)==CycleComps(idxCycleComps)-NPins)=NewComp-NPins;
    end
    C=C(2:end,:);
    if isempty(C)
        break
    end
end
EE=unique(EE,'rows');
end