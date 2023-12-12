function [A,P,basic] = pivot_correct_xpn(A,P,basic,d)
% this function makes x variables feasible
xb = basic(basic<(d+1)); % x's 
[~,pos] = ismember(xb,basic); % position of x's in basic
neg_pos = pos(A(pos,end)<0); % position of infeasible x's
for pivotRow = neg_pos
    x_pn = basic(pivotRow);
    if x_pn > d/2 % d is the number of x variables, the dimension is d/2
        pivotColumn = x_pn-d/2;
    else
        pivotColumn = x_pn+d/2;
    end
%     [A1,P1] = pivotAP(A,P,pivotRow,pivotColumn);
    A(pivotRow,:) = -A(pivotRow,:);
    P(pivotRow,:) = -P(pivotRow,:);
    basic(pivotRow) = pivotColumn;
end