function T = landmarksRelativeAngles(y)
T = zeros(1,size(y,2));
for i=2:size(y,2)
    a = y(:,i)-y(:,1);
    T(i) = a(2)/a(1);
end
end  