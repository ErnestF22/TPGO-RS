function structureMatrix_test
C=structureMatrix3();
E=eye(3);
for r=1:3
    for c=1:3
        er=E(:,r);
        ec=E(:,c);
        erc=zeros(3,1);
        idx=C(r,c);
        if idx~=0
            erc(abs(idx))=sign(idx);
        end
        if norm(erc-cross(er,ec))>0
            error('Matrix not correct')
        end
    end
end
disp('Matrix correct')
