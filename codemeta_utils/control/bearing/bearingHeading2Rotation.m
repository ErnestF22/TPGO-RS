function R=bearingHeading2Rotation(t)
if size(t,1)==2 && size(t,3)==1
    %t is really a heading vector, so convert to angle
    t=atan2(t(2,:),t(1,:));
end

R=rot_exp(eye(2),rot_hat(eye(2),t));
