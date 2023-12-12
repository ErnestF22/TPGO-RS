function y_modify = modify_bearing(L,y,idx)
T = landmarksRelativeAngles(L); %bearing between landmarks T=mxm with m landmarks
% T(i,j) is beta^j_i in the paper
ref = y(:,idx); %reference landmark
M = y(2,:)./y(1,:); %bearing from camera, it's beta_i in the paper and 
y_modify = y;
for i=1:size(y,2)
    if i~=idx% skip the one we assumed as fixed
        % find the intersection of two lines P1 and P2 such that:
        % P1: y = M(i)*x, this line passes through landmark_i
        % P2: y = T(idx,i)*x+b, this line passes through the fixed node
        b = ref(2)-T(idx,i)*ref(1);
        y_modify(:,i) = [b/(M(i)-T(idx,i));M(i)*b/(M(i)-T(idx,i))];
    end
end
end