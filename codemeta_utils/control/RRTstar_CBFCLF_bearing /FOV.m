function y_m = FOV(L,visible_set,y_visible)
% compute the bearing measurements between landmarks
T = landmarksRelativeAngles(L);
y_modify = modify_bearing(L,y_visible,visible_set(1));
y_m = zeros(2,size(L,2));
y_m(:,visible_set(1))=y_visible(:,1);
y_m(:,visible_set(2))=y_modify(:,visible_set(2));
for i=1:size(L,2)
    if ~any(i==visible_set)
        t1 = T(visible_set(1),i);
        t2 = T(visible_set(2),i);
        P = [t1, y_modify(2,1)-t1*y_modify(1,1)];
        Q = [t2, y_modify(2,2)-t2*y_modify(1,2)];
        x = (Q(2)-P(2))/(P(1)-Q(1));
        y = P(1)*x+P(2);
        y_m(:,i)=[x;y];
    end
end
end