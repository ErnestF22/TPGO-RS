function C=multikron(A,B)
szA=size(A);
szB=size(B);

NRows=szA(1)*szB(1);
NCols=szA(2)*szB(2);
idxRow=reshape(1:NRows,szB(1),szA(1));
idxCol=reshape(1:NCols,szB(2),szA(2));

C=zeros(NRows,NCols,max(size(A,3),size(B,3)));
for iRow=1:szA(1)
    for iCol=1:szA(2)
        C(idxRow(:,iRow),idxCol(:,iCol),:)=multiprod(A(iRow,iCol,:),B);
    end
end
