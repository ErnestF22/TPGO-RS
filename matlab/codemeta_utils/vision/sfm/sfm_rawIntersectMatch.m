%Returns the set of matches appearing in both sets and their indeces
%function [m,idx1,idx2]=sfm_rawIntersectMatch(m1,m2)
function [m,idx1,idx2]=sfm_rawIntersectMatch(m1,m2)
flagCollectSecondIndeces=(nargout>2);

Nm1=size(m1,2);
flagM1=false(1,Nm1);
if flagCollectSecondIndeces
    Nm2=size(m2,2);
    flagM2=false(1,Nm2);
end

for im1=1:Nm1
    id1=m1(1,im1);
    id2=m1(2,im1);
    
    if id2==sfm_rawGetMatchFromId(m2,id1)
        flagM1(im1)=true;
        if flagCollectSecondIndeces
            flagM2(id2)=true;
        end
    end
end

m=m1(:,flagM1);
if nargout>1
    idx1=find(flagM1);
    if flagCollectSecondIndeces
        idx2=find(flagM2);
    end
end

