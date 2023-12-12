function [A,basic] = pivot_correct_xpn_A(A,basic,d)
% this function makes x variables feasible
xb = basic(basic<(d+1)); % x's 
[~,pos] = ismember(xb,basic); % position of x's in basic
neg_pos = pos(A(pos,end)<0); % position of infeasible x's
for pivotRow = neg_pos
    x_pn = basic(pivotRow);
    if x_pn > d/2
        pivotColumn = x_pn-d/2;
    else
        pivotColumn = x_pn+d/2;
    end
    A(pivotRow,:) = -A(pivotRow,:);
    basic(pivotRow) = pivotColumn;
end