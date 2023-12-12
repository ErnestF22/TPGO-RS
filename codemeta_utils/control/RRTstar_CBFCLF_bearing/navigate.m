function k_path = navigate(start,tree,obs,n,st,l1,l2,l3,Landmarks,flag)
x_n=start;
u=[100;100];
e = 0.5;
e = 5; %0.5;
c=0;

k_idx  = find_nearest_k(tree,x_n,obs);
k_path = k_idx;
while norm(u)>10^-3
    c = c+1;
    x_p=x_n;
    if flag %here is just for MATLAB experiment
        L = tree(k_idx).L;
        if L == 1
            L = l1;
        elseif L==2
            L = l2;
        elseif L==3
            L = l3;
        else
            L = l4;
        end
    else %here is for Lab experiment:
        idx_landmarks = tree(k_idx).L;
        L = Landmarks(2:3,ismember(Landmarks(1,:),idx_landmarks));
        
    end
    k = tree(k_idx).K;
    Y=(L-x_p*ones(1,size(L,2)));% displacement measurements
    Y = bearing(Y,0);%notmalized the measuremnets to convert them into the bearings
%     Y = FOV(L,[1,2],Y(:,1:2)); % if you have the limited field of view
%     run this line and comment line 33, I assume landmarks 1 and 2 are
%     visible, if you have other landmarks as visible, change the index
%     accordingly
    Y = modify_bearing(L,Y,1);% if you don't have the limited field of view run this line and comment line 32
    Y = reshape(Y,[],1);% vectorized the bearing measurements
    dis = norm(x_n-tree(tree(k_idx).parent).position);% its just a threshold, you can work with the velocity
    u=k*Y;
    if dis>5%3
        u = u/norm(u);
    end
    x_n=x_p+e*u;
    if (dis<5) && ((tree(k_idx).parent)~=1)
        k_idx  = tree(k_idx).parent;
        k_path = [k_path k_idx];
    else
        if norm(u)<0.1 %0.02
            break
        else
            norm (x_n-x_p)
            if norm (x_n)<0.2
                break
            end
        end
    end
    figure(n)
    plot(x_n(1),x_n(2),st,'MarkerSize',2)
    hold on
%     figure(3)
%     plot(c,dis,st)
%     hold on
end

end
