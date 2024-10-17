function [EE,C,membership,flag_removed]=Factorgraph_removeRectangle(EE,C,membership,x,pins)
pinAddress=find(pins);
NPins=max(EE(:,1));
flag_removed=true;
while sum(C(1,:))==8
    %disp('Rectangle found');
    Csorted=sort(C(1,:));
    Va=x(:,pinAddress(Csorted(1)));
    Vb=x(:,pinAddress(Csorted(2)));
    Vc=x(:,pinAddress(Csorted(3)));
    Vd=x(:,pinAddress(Csorted(4)));
    if rank([Vb-Va,Vc-Va,Vd-Va])<3
        flag_removed=false;
        % Should I rigidify it here? probably no
        break
    end
    
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