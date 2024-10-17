function A=adjmetropolis(A)
    N=size(A,1);
    d=sum(A,2);
    for(inode=1:N)
        for(jnode=1:N)
            if(A(inode,jnode)~=0)
               A(inode,jnode)=1/(1+max(d(inode),d(jnode)));
            end
        end
    end
