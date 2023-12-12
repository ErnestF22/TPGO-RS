function landmarks = find_near_landmarks(x,y)
for i=1:size(y,2)
    dis(i) = norm(x-y(:,i));
end
[~,idx] = mink(dis,4);
landmarks = y(:,idx);
end