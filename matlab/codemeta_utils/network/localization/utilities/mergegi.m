function [t_node]=mergegi(t_node,membergi,memberRi, memberTi)
    if ~exist('membergi','var') || isempty(membergi)
        membergi='gi';
    end
    if ~exist('memberRi','var') || isempty(memberRi)
        memberRi='Ri';
    end
    if ~exist('memberTi','var') || isempty(memberTi)
        memberTi='Ti';
    end
    
    structType=testNetworkDetectStructType(t_node);

    switch structType
        case 'single'
            t_node.(membergi)=RT2G(t_node.(memberRi),t_node.(memberTi));
        case 'array'
            N=length(t_node);
            for inode=1:N
                t_node(inode).(membergi)=RT2G(t_node(inode).(memberRi),t_node(inode).(memberTi));
            end
    end
end
