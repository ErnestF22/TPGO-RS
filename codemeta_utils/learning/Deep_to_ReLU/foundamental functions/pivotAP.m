function [A,P] = pivotAP(A,P,pivotRow,pivotColumn)
ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
for i=1:(size(A,1)) %each row:
    if A(i,pivotColumn)~=0
        if i~=pivotRow
            P(i,:)=P(i,:)+ratio(i)*P(pivotRow,:);
            A(i,:)=A(i,:)+ratio(i)*A(pivotRow,:);
%             P(i,:)=P(i,:)-A(i,pivotColumn)/A(pivotRow,pivotColumn)*P(pivotRow,:);
%             A(i,:)=A(i,:)-A(i,pivotColumn)/A(pivotRow,pivotColumn)*A(pivotRow,:);

        end
    end
end
P(pivotRow,:)=P(pivotRow,:)/A(pivotRow,pivotColumn);
A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);


% ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
% ratio(pivotRow) = 1;
% PivotMatrix = eye(size(P,1));
% PivotMatrix(:,pivotRow)=ratio;
% P = PivotMatrix*P;
% A = PivotMatrix*A;
% P(pivotRow,:)=P(pivotRow,:)/A(pivotRow,pivotColumn);
% A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);

end