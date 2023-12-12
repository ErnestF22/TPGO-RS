function [basic,result,A,P]=dual_simplex2(A,basic,d)
flagFeasible = all(A(1:end-1,end)>=0);
rows = (1:size(A,1)-1);
P = eye(size(A,1));
i = 0;
result=zeros(1,size(A,2)-1);
while i<1000
    if flagFeasible
        break
    end
    neg = A(1:end-1,end)<0;
    idx_negB = rows(neg);
    pivotRow = idx_negB(1);
    if all(A(pivotRow,1:end-1)>=0)
        A = [];
        basic = [];
        result = [];
        break
    end
    flagPivotRowNegative=A(pivotRow,1:end-1)<0;
    x = zeros(1,size(A,2)-1);
    x(flagPivotRowNegative) = abs(A(pivotRow,flagPivotRowNegative));
    x(~flagPivotRowNegative) = -1;
    % exclude x+ x- in basic simultaneously
    xp_inbasic = basic(basic<=d/2);
    xn_inbasic = basic(d/2<basic & basic<=d);
    x(xp_inbasic+d/2) = -Inf;
    x(xn_inbasic-d/2) = -Inf;
    % x in basic  = -Inf
    x(basic(pivotRow)) = -Inf;
    % A==0 = -Inf
    A_piv = A(pivotRow,1:end-1);
    x(A_piv==0) = -Inf;
    [~ ,pivotCol] = max(x);
    ratio=-A(:,pivotCol)/A(pivotRow,pivotCol);
    for i=1:(size(A,1)) %each row:
        if A(i,pivotCol)~=0
            if i~=pivotRow
                A(i,:)=A(i,:)+ratio(i)*A(pivotRow,:);
                P(i,:)=P(i,:)+ratio(i)*P(pivotRow,:);
            end
        end
    end
    P(pivotRow,:)=P(pivotRow,:)/A(pivotRow,pivotCol);
    A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotCol);
    A(abs(A)<0.000001)=0;
    %update indexes in the basic set
    basic(pivotRow)=pivotCol;
    flagFeasible = all(A(1:end-1,end)>=0);
    i=i+1;
end
if ~isempty(basic)
    result(1,basic(:))=A(1:(size(A,1)-1),end);%A(1:(row-1),end) 
end
end