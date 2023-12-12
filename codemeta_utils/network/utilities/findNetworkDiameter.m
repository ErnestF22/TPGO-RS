%Computes the network diameter from the adjacency matrix A (trivial
%implementation)
function d=findNetworkDiameter(A)
    N=size(A,1);
    d=0;
    
    if(min(sum(A,2))==0)
        d=Inf;
        return
    end
    
    for(inode=1:N)
        d=max(d,findNewtorkDiameterNode(A,inode));
    end
    
end

function di=findNewtorkDiameterNode(A,inode)
    itMax=1e6;
    N=size(A,1);
    
    flagNode=zeros(N,1);
    flagNode(inode)=1;
    for(it=1:itMax)
        idx=find(flagNode~=0);
        for(jnode=1:length(idx))
            flagNode(A(idx(jnode),:)~=0)=1;
        end
        if(all(flagNode))
            di=it;
            break
        end
    end
end
