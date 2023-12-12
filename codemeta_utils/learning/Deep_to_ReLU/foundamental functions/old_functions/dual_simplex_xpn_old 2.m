function [basic,result,A,P]=dual_simplex_xpn_old(A,basic,d)
%initialization:
result=zeros(1,size(A,2)-1);
P = eye(size(A,1));
rows = (1:size(A,1)-1);

[A,P,basic] = pivot_correct_xpn(A,P,basic,d);

flagFeasible = all(A(1:end-1,end)>=0); % check initial slack feasible
flagNotFind = false; % assume we can find at least one feasible vertex
PreviousBasic = zeros(1,size(A,1)-1);
PreviousX = zeros(1,size(A,1)-1);
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
    
%     A = round(A,11);
    flagFeasible = all(A(1:end-1,end)>=0); 
end

% make sure x in basic
xb = basic(basic<(d+1)); %x's in basic
if size(xb,2)>=d/2 % at least d x's in basic
    flagXinBasic = true;
else
    flagXinBasic = false;
end
        
while (~flagXinBasic && ~flagNotFind)
    
    % LoopAllSlack
    xs = basic(basic>d); %find slack variables in the basis
    % loop all possible slack variables
    flagLoopAllSlack = true;
    for i=1:size(xs,2)
        pivotRow = find(basic==xs(i)); %row of the first slack variable
        if all(A(pivotRow,1:d)<=0) % all negative, cannot pivot
            continue
        else
            flagPivotSuccess = false;
            for j=1:d % try every x's
%                 if ismember(j,basic) % if x already in basic
%                     continue
%                 else
%                     pivotColumn = j;
%                 end
                [~,pivotColumn] = max(A(pivotRow,1:d));
                if abs(A(pivotRow,pivotColumn))<=1e-11 % cannot find a solution
                    continue
                end
                basic_new = basic;
                basic_new(pivotRow)=pivotColumn;
                % try to pivot
                [Anew,Pnew] = pivotAP(A,P,pivotRow,pivotColumn);
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
            if flagPivotSuccess
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
    idx = 1:d;
    pivotColumnSet = idx(~ismember(idx,xb_all));
    
    flagLoopAllBasic = true;
    if flagLoopAllSlack && ~isempty(pivotColumnSet)
        for i=1:size(xb,2)
            pivotRow = find(basic==xb(i)); %row of the first slack variable
            if all(A(pivotRow,pivotColumnSet)<=0) % all negative, cannot pivot
                continue
            else
                flagPivotSuccess = false;
                for j = pivotColumnSet % try every x's
                    pivotColumn = j;
                    basic_new = basic;
                    basic_new(pivotRow)=pivotColumn;
                    if abs(A(pivotRow,pivotColumn))<=1e-5 || ismember(sort(basic_new,'ascend'),PreviousBasic,'rows')% cannot find a solution
                        continue
                    end
                    % try to pivot
                    [Anew,Pnew] = pivotAP(A,P,pivotRow,pivotColumn);
                    [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
                    if all(Anew(1:end-1,end)>=0)% new solution is feasible
                        A = Anew;
                        P = Pnew;
                        basic = basic_new;
                        PreviousBasic = [PreviousBasic;sort(basic_new,'ascend')];
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
    
%     % Loop Between Basic
%     xb = basic(basic<(d+1));
%     % find corresponding x+,x- except themselves
%     xb_all = [xb+d/2,xb-d/2];
%     xb_all(xb_all>d) = [];
%     xb_all(xb_all<1) = [];
%     % columns that can be pivoted to (x+ to x-, or x- to x+)
% %     pivotColumnSet = xb_all;
%     pivotColumnSet = [pivotColumnSet,xb_all];
%     
%     flagLoopBetweenBasic = true;
%     if flagLoopAllSlack && flagLoopAllBasic %&& 0 %%%%%%%%%%
%         for i=1:size(xb,2)
%             pivotRow = find(basic==xb(i)); %row of the first slack variable
%             flagPivotSuccess = false;
%             for j = pivotColumnSet % try every x's
%                 pivotColumn = j;
%                 if abs(A(pivotRow,pivotColumn))<=1e-11 % cannot find a solution
%                     continue
%                 end
%                 % try to pivot
%                 [FlagPivotFeasible,Anew,Pnew,basic_new] = pivotfeasible(A,P,pivotRow,pivotColumn,basic);
%                 if FlagPivotFeasible && ~ismember(sort(basic_new,'ascend'),PreviousX,'rows')% new solution is feasible
%                     A = Anew;
%                     P = Pnew;
%                     basic = basic_new;
%                     PreviousX = [PreviousX;sort(basic_new,'ascend')];
%                     flagLoopBetweenBasic = false; % find a pivotable slack variable, no need to loop all
%                     flagPivotSuccess = true;
%                     break
%                 else
%                     continue
%                 end
%             end
%             if flagPivotSuccess
%                 break
%             end
%         end
%     end
    
    flagLoopBetweenBasic = true;
    if flagLoopAllSlack  && flagLoopAllBasic && flagLoopBetweenBasic % loop all, but cannot find
        flagNotFind = true;
        break
    end
    
    xb = basic(basic<(d+1)); %x's in basic (2*d+1)
    if size(xb,2)>=d/2 %at least d x's in basic
        flagXinBasic = true;
        xb_check = [xb+d/2,xb-d/2];
        xb_check(xb_check>d) = [];
        xb_check(xb_check<1) = [];
        [val,~]=intersect(xb,xb_check);
        if ~isempty(val)
            disp('repeat!')
        end
    end
end

% new basic
if flagNotFind % if do not find a solution
    A = [];
    basic = [];
    result = [];
%     disp('Vertex Not Find');
else
    result(1,basic(:))=A(1:(size(A,1)-1),end);
end
end