function Q = pivot_flip_matrix(basic,A,P)
Q(1).basic = basic;
Q(1).A = A;
Q(1).P = P;
d = size(A,2)-size(A,1);
for xs = basic % for every variable in basic
    Arf = A;
    Pnew = P;
    if xs<=d % do not pivot for x variables
        continue
    end
    pivotRow = find(basic==xs); 
    % find pivot column
    for i=(d+1):size(A,2)-1 %for every column in A except x variables
        if any(basic==i) % if already in basic
            continue
        else
            % find new base
            pivotColumn = i;
            basic_new = basic;
            basic_new(basic_new == xs)=i;
            if Arf(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                continue
            else % pivot
                [Anew,Pnew] = pivotAP(Arf,Pnew,pivotRow,pivotColumn);
                [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
                if any(Anew(1:end-1,end)<0) % not a feasible solution
                    continue
                else
                    Q(end+1).basic = basic_new;
                    Q(end).A = Anew;
                    Q(end).P = Pnew;
                end
            end
        end
    end
end
end