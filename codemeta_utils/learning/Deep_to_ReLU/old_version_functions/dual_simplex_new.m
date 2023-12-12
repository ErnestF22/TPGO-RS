function [basic,result,A,P]=dual_simplex_new(A,basic,d)
%initialization:
columns=size(A,2);
result=zeros(1,size(A,2)-1);
flagBasicSolutionFeasible=false;
flagVer = false;
P = eye(size(A,1));
x=zeros(1,columns);
x(1,columns)=Inf;
% while flagBasicSolutionInfeasible
while ~(flagVer && flagBasicSolutionFeasible)
    if min(A(1:(size(A,1)-1),end))<0 %if there is infeasible solution
        [~ ,pivotRow]=min(A(1:(size(A,1)-1),end));    %A(1:row-1,end) minimum of last column
        %find pivot row
        flagAPivotRowNegative=A(pivotRow,1:end-1)<0;
        x(flagAPivotRowNegative)=A(end,flagAPivotRowNegative)./abs(A(pivotRow,flagAPivotRowNegative));
        x(~flagAPivotRowNegative)=Inf;
        [~ ,pivotColumn]=min(x);
    else
        xs = basic(basic>(2*d)); %find slack variables in the basis
        pivotRow = find(basic==xs(1)); %row of the first slack variable
        [~,pivotColumn] = max(A(pivotRow,1:(2*d)));
    end

    if abs(A(pivotRow,pivotColumn)) <= 1e-10
        A = [];
        result=[];
        basic=[];
        break
    end
    %perform pivot (transform A such that pivot column becomes
    %[0;0;...;0;1;0;...;0]
%     for i=1:(size(A,1)) %each row:
%         if A(i,pivotColumn)~=0
%             if i~=pivotRow
%                 ratio=-A(i,pivotColumn)/A(pivotRow,pivotColumn);
%                 A(i,:)=A(i,:)+ratio*A(pivotRow,:);
%             else
%                 A(pivotRow,:)=A(pivotRow,:)/A(pivotRow,pivotColumn);
%             end
%         end
%     end
    
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
    
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;

    %if the basic solution (last column of A) is feasible, we will stop
    flagBasicSolutionFeasible=~any(A(1:end-1,end)<0);
    
    xb = basic(basic<(2*d+1)); %x's in basic (2*d+1)
    if size(xb,2)>=d  %at least two x's in basic (d)
        flagVer = true;
    end
    
%     A = round(A,11);
%     P = round(P,10);

end
% new basic
if ~isempty(basic)
    result(1,basic(:))=A(1:(size(A,1)-1),end);    %A(1:(row-1),end) 
end

end