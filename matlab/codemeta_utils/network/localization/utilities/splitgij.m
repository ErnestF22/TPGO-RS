function [t_node]=splitgij(t_node,membergij,memberRij, memberTij)
    if ~exist('membergij','var') || isempty(membergij)
        membergij='gij';
    end
    if ~exist('memberRij','var') || isempty(memberRij)
        memberRij='Rij';
    end
    if ~exist('memberTij','var') || isempty(memberTij)
        memberTij='Tij';
    end
    
    structType=testNetworkDetectStructType(t_node);

    switch structType
        case 'single'
            [t_node.(memberRij),t_node.(memberTij)]=G2RT(t_node.(membergij));
        case 'array'
            N=length(t_node);
            for inode=1:N
                for jnode=1:N
                    [t_node(inode).(memberRij)(:,:,jnode),t_node(inode).(memberTij)(:,jnode)]=...
                        G2RT(t_node(inode).(membergij)(:,:,jnode));
                end
            end
    end
end

