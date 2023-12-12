function [basic,result,A,P]=dual_simplex1(A,basic,d)

%initialization:
columns=size(A,2);
result=zeros(1,size(A,2)-1);
flagBasicSolutionFeasible=false;
flagVer = false;
P = eye(size(A,1));

% while flagBasicSolutionInfeasible
while ~( flagVer && flagBasicSolutionFeasible)
    neg = A(1:(size(A,1)-1),end)<0;
    neg_b = basic(neg);
    Slack_Infeasible = any(neg_b>d);
    if Slack_Infeasible %if there is infeasible solution
        pivotidx = neg_b(find(neg_b>d,1));
        pivotRow = find(basic==pivotidx);
        %find pivot row
        flagAPivotRowNegative=A(pivotRow,1:end-1)<0;
        x=zeros(1,columns);
        x(1,columns)=Inf;
        x(flagAPivotRowNegative)=A(end,flagAPivotRowNegative)./abs(A(pivotRow,flagAPivotRowNegative));
        x(~flagAPivotRowNegative)=Inf;
        [~ ,pivotColumn]=min(x);
    else
        xs = basic(basic>d); %find slack variables in the basis
        pivotRow = find(basic==xs(1)); %row of the first slack variable
        [~,pivotColumn] = max(abs(A(pivotRow,1:d)));
    end
    
    if abs(A(pivotRow,pivotColumn)) <= 1e-10
        A = [];
        result=[];
        basic=[];
        break
    end
    %perform pivot (transform A such that pivot column becomes
    %[0;0;...;0;1;0;...;0]
    % include pivoting for P
    ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
    for i=1:(size(A,1)) %each row:
        if A(i,pivotColumn)~=0
            if i~=pivotRow
                A(i,:)=A(i,:)+ratio(i)*A(pivotRow,:);
                P(i,:)=P(i,:)+ratio(i)*P(pivotRow,:);
            end
        end
    end
    P(pivotRow,:)=P(pivotRow,:)/A(pivotRow,pivotColumn);
    A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);
    A(abs(A)<0.000001)=0;
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;

    %if the basic solution (last column of A) is feasible, we will stop
    neg = A(1:(size(A,1)-1),end)<0;
    neg_b = basic(neg);
    flagBasicSolutionFeasible = ~any(neg_b>d);
    flagBasicSolutionFeasible = all(A(1:(size(A,1)-1),end)>=0);
    xb = basic(basic<(d+1)); %x's in basic (2*d+1)
    if size(xb,2)>=d %at least d x's in basic (d)
        flagVer = true;
    end
end
% new basic
if ~isempty(basic)
    result(1,basic(:))=A(1:(size(A,1)-1),end);%A(1:(row-1),end) 
end

end

%         D = (1:d);
%         xx = basic(basic<=d);
%         xa = D(~ismember(D,xx));
%         C = A(1:end-1,xa(1))<0;
%         Cb = basic(C);