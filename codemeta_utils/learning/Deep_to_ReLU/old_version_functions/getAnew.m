function [A,basic] = getAnew(A,act,basic)
d = size(A,2)-size(A,1);
rows = (1:size(A,1)-1);
[A,basic] = pivot_correct_xpn_A(A,basic,d);
flagFeasible = all(A(1:end-1,end)>=0); % check initial slack feasible
flagNotFind = false;

while (~flagFeasible && ~flagNotFind) % pivot toward a feasible vertex
    neg = A(1:end-1,end)<0;
    idx_negB = rows(neg);
    pivotRow = idx_negB(1);
    if all(A(pivotRow,1:end-1)>=0) % all positive, cannot pivot -> cannot find a solution
        flagNotFind = true;
        break
    end
    flagPivotRowNegative=A(pivotRow,1:end-1)<0;
    x = zeros(1,size(A,2)-1);
    x(flagPivotRowNegative) = abs(A(pivotRow,flagPivotRowNegative));
    x(~flagPivotRowNegative) = -Inf;
    [~ ,pivotColumn] = max(x);% maximum of the abs(negative values)
    
    if abs(A(pivotRow,pivotColumn))<=1e-11 % cannot find a solution
        flagNotFind = true;
        break
    end
    
    % include pivoting for P
    A = pivotA(A,pivotRow,pivotColumn);
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;
    
    [A,basic] = pivot_correct_xpn_A(A,basic,d);
    
    flagFeasible = all(A(1:end-1,end)>=0); 
end

xb = basic(basic<=d);
if flagNotFind
    A = [];
    basic = [];
elseif size(xb,2)==d/2
    
    for pivotRow = act
        if ismember(pivotRow,basic)
            pivotRow = find(basic==pivotRow);
        else
            continue
        end
        xb = basic(basic<=d);
        for pivotColumn = 1:d/2
            if ismember(pivotColumn,xb)
                continue
            else
                basic_new = basic;
                A_new = A;
                A_new = pivotA(A_new,pivotRow,pivotColumn);
                basic_new(pivotRow) = pivotColumn;
                if any(A_new(1:end-1,end)<0)
                    [A_new,basic_new] = pivot_correct_xpn_A(A_new,basic_new,d);
                end
                if any(A_new(1:end-1,end)<0)
                    continue
                else
                    basic = basic_new;
                    A = A_new;
                    break
                end
            end
        end 
    end
end