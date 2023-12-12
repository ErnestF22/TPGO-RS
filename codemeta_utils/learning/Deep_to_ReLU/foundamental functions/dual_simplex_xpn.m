function [basic,A,P]=dual_simplex_xpn(A,basic)
%initialization:
result=zeros(1,size(A,2)-1);
d = size(A,2)-size(A,1);
P = eye(size(A,1));
rows = (1:size(A,1)-1);

[A,P,basic] = pivot_correct_xpn(A,P,basic,d);

flagFeasible = all(A(1:end-1,end)>=0); % check initial slack feasible
flagNotFind = false; % assume we can find at least one feasible vertex
PreviousBasic = zeros(1,size(A,1)-1); % save previously checked basic

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
    [A,P] = pivotAP(A,P,pivotRow,pivotColumn);
    %update indexes in the basic set
    basic(pivotRow)=pivotColumn;
    
    [A,P,basic] = pivot_correct_xpn(A,P,basic,d);
    
    flagFeasible = all(A(1:end-1,end)>=0); 
end

% make sure x in basic
xb = basic(basic<(d+1)); %x's in basic
if size(xb,2)>=d/2 % at least d x's in basic
    flagXinBasic = true;
else
    flagXinBasic = false;
end
% initialize all x's
x_idx = 1:d;
        
while (~flagXinBasic && ~flagNotFind)
    % LoopAllSlack
    xs = basic(basic>d); %find slack variables in the basis
    % loop all possible slack variables
    flagLoopAllSlack = true;
    for i=1:size(xs,2)
        pivotRow = find(basic==xs(i)); %row of the first slack variable
        if all(A(pivotRow,1:d)<=0) % all x's negative, cannot pivot
            continue
        else
            flagPivotSuccess = false;
            pivotColumnSet = find(A(pivotRow,1:d)>0); % pivotable column
            sizeColumn = size(pivotColumnSet,2);
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
                    flagLoopAllSlack = false; % find a pivotable slack variable, no need to loop all
                    flagPivotSuccess = true;
                    break
                else
                    continue
                end
            end
            if flagPivotSuccess % if find feasible
                break
            end
        end
    end
    
    % LoopAllBasic
    xb = basic(basic<(d+1));
    % find corresponding x+,x-
    xb_all = [xb,xb+d/2,xb-d/2];
    xb_all(xb_all>d) = [];
    xb_all(xb_all<1) = [];
    % columns that can be pivoted to
    pivotColumnSet = x_idx(~ismember(x_idx,xb_all));
    % loop all possible x's variables
    flagLoopAllBasic = true;
    if flagLoopAllSlack && ~isempty(pivotColumnSet)
        for i=1:size(xb,2)
            pivotRow = find(basic==xb(i)); %row of the first x variable
            if all(A(pivotRow,pivotColumnSet)<=0) % all negative, cannot pivot
                continue
            else
                flagPivotSuccess = false;
                sizeColumn = size(pivotColumnSet,2);
                for j=1:sizeColumn % try every x's
                    [~,pivotColumnPosition] = max(A(pivotRow,pivotColumnSet));
                    pivotColumn = pivotColumnSet(pivotColumnPosition);
                    pivotColumnSet(pivotColumnSet == pivotColumn) = [];
                    basic_new = basic;
                    basic_new(pivotRow)=pivotColumn;
                    if abs(A(pivotRow,pivotColumn))<=1e-5 || sum(sum(PreviousBasic-sort(basic_new),2)==0)~=0% cannot find a solution
                        continue
                    end
                    % try to pivot
                    [Anew,Pnew] = pivotAP(A,P,pivotRow,pivotColumn);
                    [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
                    if all(Anew(1:end-1,end)>=0)% new solution is feasible
                        A = Anew;
                        P = Pnew;
                        basic = basic_new;
                        PreviousBasic = [PreviousBasic;sort(basic_new)]; % store checked basic to avoid falling into loops
                        flagLoopAllBasic = false; % find a pivotable slack variable, no need to loop all
                        flagPivotSuccess = true;
                        break
                    else
                        continue
                    end
                end
                if flagPivotSuccess
                    break
                end
            end
        end
    end
    
    if flagLoopAllSlack  && flagLoopAllBasic % loop all, but cannot find
        flagNotFind = true;
        break
    end
    
    xb = basic(basic<(d+1)); %x's in basic (2*d+1)
    if size(xb,2)>=d/2 %at least d x's in basic
        flagXinBasic = true;
    end
end

% new basic
if flagNotFind % if do not find a solution
    A = [];
    basic = [];
    result = [];
else
    result(1,basic(:))=A(1:(size(A,1)-1),end);
end
end