% this function computes the angle between two vector:
% vector vec_np: vectro from node n to node p
% vector vec_obs: vectro from node n to obstacle point
% This function is a part of nearest_angle_Obstacle.m function
function [ang_l,ang_r] = compute_angle_obs(n,p,obs,d)
ang_l = [];
ang_r = [];
vec_np = n-p;
left = [];
right = [];
right_index = [];
left_index =[];
for i=1:size(d,2)
    if d(i)~=Inf
        vec_obs = n - obs(:,i);
        z = [0 0 1]*cross([vec_np;0],[vec_obs;0]);
        t = acos(dot(vec_np,vec_obs)/(norm(vec_np)*norm(vec_obs)));
        if z>=0
            left = [left t];
            left_index = [left_index i];
        else
            right = [right t];
            right_index = [right_index i];
        end
    end
end
if ~isempty(right)
[~,idx] = min(right);
ang_r = right_index(idx);
end
if ~isempty(left)
[~,idx] = min(left);
ang_l = left_index(idx);
end
end