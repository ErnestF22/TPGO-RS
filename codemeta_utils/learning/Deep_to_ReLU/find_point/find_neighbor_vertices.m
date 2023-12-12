function [Q,ACT] = find_neighbor_vertices(Q,min_idx,ACT,x,Ab_set,plane,flagCheck)
% find the neighbor vertices of one vertex
A = Q(min_idx).A;
basic = Q(min_idx).basic;
act = Q(min_idx).act;
z = Q(min_idx).z;
P = Q(min_idx).P;
size_var = size(A,2)-1; % number of all variables

d = size(A,2)-size(A,1);
zflip_set = unique(act)-d; % z index to flip from active constraints
zflip_set = [0,zflip_set]; % 0 index means not flipping(inside the region)

A_add = plane.A;
b_add = plane.b;

for zflip = zflip_set
    if zflip == 0 % find the neighboring vertices inside current region
        znew = z;
        Arf = A;
        Prf = P;
    else
        znew = z;
        znew(zflip) = ~znew(zflip); % z after flipping
        [Arf,Prf] = find_tableau_after_flip_P(z,zflip,P,A,Ab_set); % tableau and P after flipping
        if flagCheck
            [Q1,ACTset] = checkA(z,zflip,Ab_set,x); % use dual simplex to find vertices
            [~,col] = ismember(act,ACTset,'rows'); % find the position of this vertex in Q
            sumall = abs(sum((Q1(col).A-Arf),'all')); % the differene between tableau from find_vertices and P matrix
            % test if the final matrix match
            if sumall<=1e-1
                match = ['match! difference is ',num2str(sumall)];
                disp(match);
            else
                notmatch = ['not match! difference is ',num2str(sumall)];
                disp(notmatch);
            end
            flagCheck = false;
        end
    end
    for xs = basic % for every variable in basic
        if xs<=d % do not pivot for x variables
            continue
        end
        pivotRow = find(basic==xs); 
        % find pivot column
        for i = act % for every active variables
            % find new base
            pivotColumn = i;
            basic_new = basic;
            basic_new(basic_new == xs)=i;
            a_new = act;
            a_new(a_new == pivotColumn) = xs;
            a_new_sorted = sort(a_new);
            if Arf(pivotRow,pivotColumn)<= 1e-5 % if the pivot value <=0, decided by A
                continue
            else % pivot
                if sum(sum(abs(ACT-a_new_sorted),2)==0)==0 % if ~ismember(sort(a_new),ACT,'rows')
                    [Anew,Pnew] = pivotAP(Arf,Prf,pivotRow,pivotColumn);
                    [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
                    if any(Anew(1:end-1,end)<0) % not a feasible solution
                        continue
                    else
                        result = zeros(1,size_var);
                        result(1,basic_new(:)) = Anew(1:end-1,end);
                        v = result(1:d/2)-result(d/2+1:d);
                        ACT = [ACT;a_new_sorted];
                        Q(end+1).act = a_new_sorted;
                        Q(end).basic = basic_new;
                        Q(end).v = v;
                        Q(end).A = Anew;
                        Q(end).z = znew;
                        h_input = get_h_input(x,v,Anew,a_new_sorted);
                        h_plane = abs(A_add*v'-b_add)/norm(A_add); % distance from the target plane
                        Q(end).h_input = h_input;
                        Q(end).h_plane = h_plane;
                        Q(end).P = Pnew;
                    end
                end
            end
        end
    end
end
end