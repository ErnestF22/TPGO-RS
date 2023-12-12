function normalized_y=bearing(y,e)
normalized_y = zeros(size(y,1),size(y,2));

for i=1:size(y,2)
    normalized_y(:,i) = y(:,i)/norm(y(:,i));
%     cos_t = dot(e,y(:,i))/(norm(e)*norm(y(:,i)));
%     norm_y = 1/(cos_t*norm(e));
%     N = norm(y(:,i))/norm_y;
%     normalized_y(:,i) = y(:,i)/N;
end
end