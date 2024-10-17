%function [t_node]=splitgi(t_node,membergi,memberRi, memberxi)
%Splits the poses in t_node.(membergi) into t_node.(memberRi) and
%t_node.(memberxi). If the member names are omitted or empty, these are the
%default values: membergi='gi', memberRi='Ri' and memberxi='xi'.
function [t_node]=splitgi(t_node,membergi,memberRi, memberTi)
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
            [t_node.(memberRi),t_node.(memberTi)]=G2RT(t_node.(membergi));
        case 'array'
            N=length(t_node);
            for inode=1:N
                [t_node(inode).(memberRi),t_node(inode).(memberTi)]=...
                    G2RT(t_node(inode).(membergi));
            end
    end
    
end
