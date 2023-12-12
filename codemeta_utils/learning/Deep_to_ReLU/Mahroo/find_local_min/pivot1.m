function [flagFind,A,Q,basic] = pivot1(A,pivotRow,pivotColumn,basic,d,Q)

for i=1:(size(A,1)) %each row:
        if A(i,pivotColumn)~=0
            if i~=pivotRow
                ratio=-A(i,pivotColumn)/A(pivotRow,pivotColumn);
                A(i,:)=A(i,:)+ratio*A(pivotRow,:);
            else
                a=1;
                A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);
            end
        end
end

basic(pivotRow)=pivotColumn;
    neg = find(A(1:(size(A,1)-1),end)<0);
    neg_b = basic(neg);
    flagBasicSolutionFeasible = ~any(neg_b>d);

     flagBasicSolutionFeasible= 0;
end
xb = basic(basic<(d+1)); %x's in basic (d+1)
if size(xb,2)>=d %at least d x's in basic (d)
    flagVer = true;
else
    flagVer = false;
end

flagExist = false;

columns=size(A,2);

% while flagBasicSolutionInfeasible
while ~(flagVer && flagBasicSolutionFeasible)
    neg = find(A(1:(size(A,1)-1),end)<0);
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
        basic=[];
        break
    end
    
    %perform pivot
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
    
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;
    base = sort(basic,'ascend');
    if sum(sum(Q-base,2)==0)% if already in Q
        flagExist = true;
        break
    else
        Q = [Q;base];
    end

    %if the basic solution (last column of A) is feasible, we will stop
    neg = find(A(1:(size(A,1)-1),end)<0);
    neg_b = basic(neg);
    flagBasicSolutionFeasible = ~any(neg_b>d);
    
    xb = basic(basic<(d+1)); %x's in basic (d+1)
    if size(xb,2)>=d %at least d x's in basic (d)
        flagVer = true;
    end
    
    A = round(A,11);
    
end

% new basic
if ~isempty(basic) && ~flagExist
    flagFind = true;
else
    flagFind = false;
end

end