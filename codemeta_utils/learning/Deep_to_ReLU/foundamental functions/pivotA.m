function A = pivotA(A,pivotRow,pivotColumn)
% A4 = A;
% for i=1:(size(A,1)) %each row:
%     if A(i,pivotColumn)~=0
%         if i~=pivotRow
%             A(i,:)=A(i,:)-A(i,pivotColumn)/A(pivotRow,pivotColumn)*A(pivotRow,:);
%         end
%     end
% end
% A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);

% ratio=A4(:,pivotColumn)/A4(pivotRow,pivotColumn);
% ratio(pivotRow) = 1-1/A4(pivotRow,pivotColumn);
% A1 = diag(ratio);
% A2 = repmat(A4(pivotRow,:),size(A,1),1);
% A = A4-A1*A2;

ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
ratio(pivotRow) = 1/A(pivotRow,pivotColumn)-1;
Arow = A(pivotRow,:);
for i=1:(size(A,1)) %each row:
    if A(i,pivotColumn)~=0
        A(i,:)=A(i,:)+ratio(i)*Arow;
    end
end

end