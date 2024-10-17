function [IndexVec,flagTerminate]=combnkIndexGen(n,k,v)
if k>n
    error('Wrong Inputs!');
end
if isempty(v)
    IndexVec=1:k;
    flagTerminate=false;
    return;
end
flagTerminate=false;
index=k;
upperlimit=n;

while v(index)==upperlimit
    index=index-1;
    upperlimit=upperlimit-1;
    if index==0
        IndexVec=1:k;
        flagTerminate=true;
        return;
    end
end
v(index)=v(index)+1;
for idx=index+1:k
    v(idx)=v(index)+idx-index;
end
IndexVec=v;

end