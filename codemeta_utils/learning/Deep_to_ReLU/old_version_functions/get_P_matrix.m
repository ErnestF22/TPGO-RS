function P = get_P_matrix(A,act)
d = size(A,2)-size(A,1);
P = eye(size(A,1));
rows = (1:size(A,1)-1);

s = size(A,1); % slack variables
var = d+s; % all variables
basic = (var-s+1):var;

[A,P,basic] = pivot_correct_xpn(A,P,basic,d);

flagFeasible = all(A(1:end-1,end)>=0); % check initial slack feasible

while ~flagFeasible % pivot toward a feasible vertex
    neg = A(1:end-1,end)<0;
    idx_negB = rows(neg);
    pivotRow = idx_negB(1);
    flagPivotRowNegative=A(pivotRow,1:end-1)<0;
    x = zeros(1,size(A,2)-1);
    x(flagPivotRowNegative) = abs(A(pivotRow,flagPivotRowNegative));
    x(~flagPivotRowNegative) = -Inf;
    [~ ,pivotColumn] = max(x);% maximum of the abs(negative values)
    
    % include pivoting for P
    [A,P] = pivotAP(A,P,pivotRow,pivotColumn);
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;
    
    [A,P,basic] = pivot_correct_xpn(A,P,basic,d);
    
    flagFeasible = all(A(1:end-1,end)>=0); 
end

for a = act
    pivotRow = find(basic==a);
    if isempty(pivotRow)
        continue
    end
    pivotColumnSet = find(A(pivotRow,1:d)>0); % pivotable column
    sizeColumn = size(pivotColumnSet,2);
    pivotsuccess = false;
    for j=1:sizeColumn % try every x's
        [~,pivotColumnPosition] = max(A(pivotRow,pivotColumnSet));
        pivotColumn = pivotColumnSet(pivotColumnPosition);
        pivotColumnSet(pivotColumnSet == pivotColumn) = [];
        if abs(A(pivotRow,pivotColumn))<=1e-11 % cannot find a solution
            continue
        end
        % try to pivot
        [Anew,Pnew] = pivotAP(A,P,pivotRow,pivotColumn);
        basic_new = basic;
        basic_new(pivotRow)=pivotColumn;
        [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
        if all(Anew(1:end-1,end)>=0) % new solution is feasible
            A = Anew;
            P = Pnew;
            basic = basic_new;
            pivotsuccess = true;
            break
        else
            continue
        end
    end
    if ~pivotsuccess
        disp('cannot pivot!');
    end
end
end