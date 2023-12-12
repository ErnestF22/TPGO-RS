function v = bearing_modification(v,y,x)
T = landmarksRelativeAngles(y);
for i=2:size(y,2)
    l_y = lineEqu(y(:,1),T(i));
    l_x = lineEqu(y(:,1),x);
    p  = intersectionTwoLines(l_y,l_x);
    v(:,i) = v(:,i)*norm(x-p);
end
end