% In this function compute the corresponding skew symmetric matrix of a
% given vector
function skewMatrix = hat(vec)
skewMatrix = [0 -vec(3) vec(2); vec(3) 0 -vec(1); -vec(2) vec(1) 0];
end