function A = pivot(A,pivotRow,pivotColumn)
for i=1:(size(A,1)) %each row:
        if A(i,pivotColumn)~=0
            if i~=pivotRow
                ratio=-A(i,pivotColumn)/A(pivotRow,pivotColumn);
                A(i,:)=A(i,:)+ratio*A(pivotRow,:);
            else
                A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);
            end
        end
end