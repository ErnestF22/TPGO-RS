function [flagFind,A,Q,basic] = pivot1(A,pivotRow,pivotColumn,basic,d,Q)
flagExist = false;
columns=size(A,2);
x=zeros(1,columns);
x(1,columns)=Inf;

ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
for i=1:(size(A,1)) %each row:
    if A(i,pivotColumn)~=0
        if i~=pivotRow
            A(i,:)=A(i,:)+ratio(i)*A(pivotRow,:);
        end
    end
end
A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);

basic(pivotRow)=pivotColumn;

neg_b = basic(A(1:(size(A,1)-1),end)<0);
flagBasicSolutionFeasible = ~any(neg_b>d);

xb = basic(basic<(d+1)); % x's in basic (d+1)
if size(xb,2)>=d % at least d x's in basic (d)
    flagVer = true;
else
    flagVer = false;
end

% if this is not a basic feasible vertex, continue pivoting
while ~(flagVer && flagBasicSolutionFeasible)
    if ~flagBasicSolutionFeasible %if there is infeasible solution
        pivotidx = neg_b(find(neg_b>d,1));
        pivotRow = find(basic==pivotidx);
        %find pivot row
        flagAPivotRowNegative=A(pivotRow,1:end-1)<0;
        x(flagAPivotRowNegative)=A(end,flagAPivotRowNegative)./abs(A(pivotRow,flagAPivotRowNegative));
        x(~flagAPivotRowNegative)=Inf;
        if any(x~=Inf) % if x has negative
            [~ ,pivotColumn]=min(x);
        else % if x all positive, choose basic variable to pivot
            x = -A(pivotRow,1:d);
            [~ ,pivotColumn]=min(x);
        end
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
    ratio=-A(:,pivotColumn)/A(pivotRow,pivotColumn);
    for i=1:(size(A,1)) %each row:
        if A(i,pivotColumn)~=0
            if i~=pivotRow
                A(i,:)=A(i,:)+ratio(i)*A(pivotRow,:);
            end
        end
    end
    A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);
    
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