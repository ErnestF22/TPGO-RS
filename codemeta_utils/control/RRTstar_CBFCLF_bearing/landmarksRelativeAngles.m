function T = landmarksRelativeAngles(y)
T = zeros(size(y,2),size(y,2));
for i=1:size(y,2)
    for j=1:size(y,2)
        a = y(:,j)-y(:,i);
        T(i,j) = a(2)/a(1);
    end
end
end