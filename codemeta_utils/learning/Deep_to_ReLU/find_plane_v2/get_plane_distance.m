function h_plane_min = get_plane_distance(x,z,vset,act,As,Ab_set)
d = (size(As,2)-size(As,1));
act_idx = act-d;
A1 = As(act_idx,1:d/2);
b1 = As(act_idx,end);
% this point needs to satisfy two conditions:
% 1 point on the plane
% 2 point inside region
t = (b1-A1*x)/(-sum(A1.^2));
point_v = x'-t*A1;
norm_len = abs(t*norm(A1));
r1 = A1*point_v'-b1;
if ~isequal(round(r1,5),0)
    disp('not on plane!');
end
[~,zp1] = forward(point_v',size(Ab_set,2),Ab_set); % z when input = point_v1
zp2 = zp1;
zp2(act_idx)=~zp2(act_idx); % since point_v on edge, there are two possibilities of z
if sum(abs(zp1-z))==0 || sum(abs(zp2-z))==0 % if point_v1 is inside current region
    h_plane = norm_len;
else
    h_plane = Inf;
end
% compare distance from plane and vertices
h_point_set = vecnorm(vset-x',2,2);
h_point = min(h_point_set);
if h_plane<h_point
    h_plane_min = h_plane;
else
    h_plane_min = h_point;
end
end